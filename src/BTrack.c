//=======================================================================
/** @file BTrack.c
 *  @brief BTrack - a real-time beat tracker
 *  @author Zach Banks, Adam Stark
 *  @copyright Copyright (C) 2015 Zach Banks
 *  @copyright Copyright (C) 2008-2014  Queen Mary University of London
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//=======================================================================

#include "BTrack.h"
#include "samplerate.h"
#include "common.h"


static void resampleOnsetDetectionFunction(struct btrack * bt);
static void calculateTempo(struct btrack * bt);
static void adaptiveThreshold(double * x, int N, double * aux_buffer);
static void calculateOutputOfCombFilterBank(struct btrack * bt);
static void calculateBalancedACF(double * onsetDetectionFunction, double * acf);
static double calculateMeanOfArray(double * array, int startIndex, int endIndex);
static void normaliseArray(double * array, int N);
static void updateCumulativeScore(struct btrack * bt, double odfSample);
static void predictBeat(struct btrack * bt);

// Ex. hop_size = 512; frame_size = 1024
int btrack_init(struct btrack * bt, int hop_size, int frame_size, int sample_rate){
    int rc = 0;
    bt->hopSize = hop_size;
    bt->frameSize = frame_size;
    bt->invSampleRate = 1.0 / (double) sample_rate;

    rc = odf_init(&bt->odf, hop_size, frame_size, ComplexSpectralDifferenceSq, HanningWindow); 
    BTRACK_ASSERT(rc == 0);

    double rayparam = 43;
	double pi = M_PI;
	
	// initialise parameters
	bt->tightness = 5;
	bt->alpha = 0.9;
	bt->tempo = 120;
	bt->estimatedTempo = 120.0;
	bt->tempoToLagFactor = 60.* (double) sample_rate / (double) hop_size;
	
	bt->m0 = 10;
	bt->beatCounter = -1;
	
	bt->beatDueInFrame = false;

	// create rayleigh weighting vector
	for (int n = 0;n < 128;n++) {
		bt->weightingVector[n] = ((double) n / pow(rayparam,2)) * exp((-1*pow((double)-n,2)) / (2*pow(rayparam,2)));
	}
	
	// initialise prev_delta
	for (int i = 0;i < 41;i++) {
		bt->prevDelta[i] = 1;
	}
	
	double t_mu = 41/2;
	double m_sig;
	double x;
	// create tempo transition matrix
	m_sig = 41/8;
	for (int i = 0;i < 41;i++) {
		for (int j = 0;j < 41;j++) {
			x = j+1;
			t_mu = i+1;
			bt->tempoTransitionMatrix[i][j] = (1 / (m_sig * sqrt(2*pi))) * exp( (-1*pow((x-t_mu),2)) / (2*pow(m_sig,2)) );
		}
	}
	
	// tempo is not fixed
	bt->tempoFixed = false;
    
    // initialise latest cumulative score value
    // in case it is requested before any processing takes place
    bt->latestCumulativeScoreValue = 0;
    bt->latestODF = 0;

    btrack_set_hop_size(bt, hop_size);
    return 0;
}

void btrack_del(struct btrack * bt){
    free(bt->onsetDF);
    free(bt->cumulativeScore);
    odf_del(&bt->odf);
}

int btrack_beat_due_in_current_frame(const struct btrack * bt){
    return bt->beatDueInFrame;
}

double btrack_get_bpm(const struct btrack * bt){
    return bt->estimatedTempo;
}

double btrack_get_latest_score(const struct btrack * bt){
    return bt->latestCumulativeScoreValue;
}

double btrack_get_latest_odf(const struct btrack * bt){
    return bt->latestODF;
}

double btrack_get_latest_confidence(const struct btrack * bt){
    return bt->latestConfidence;
}

int btrack_get_frames_until_beat(const struct btrack * bt) {
    return bt->beatCounter;
}

double btrack_get_time_until_beat(const struct btrack * bt) {
    double n_frames = btrack_get_frames_until_beat(bt);
    return n_frames * bt->hopSize * bt->invSampleRate;
}

void btrack_process_audio_frame(struct btrack * bt, const btrack_chunk_t * frame){
    double sample = odf_process_frame(&bt->odf, frame);
    btrack_process_odf_sample(bt, sample);
}

void btrack_process_fft_frame(struct btrack * bt, const btrack_chunk_t * frame){
    double sample = odf_process_fft_frame(&bt->odf, frame);
    btrack_process_odf_sample(bt, sample);
}

void btrack_process_odf_sample(struct btrack * bt, double odf_sample){
    // we need to ensure that the onset
    // detection function sample is positive
    // add a tiny constant to the sample to stop it from ever going
    // to zero. this is to avoid problems further down the line
    odf_sample = fabs(odf_sample) + 0.0001;
    
	bt->m0--;
	bt->beatCounter--;
	bt->beatDueInFrame = false;
	
	// move all samples back one step
	for (int i=0;i < (bt->onsetDFBufferSize-1);i++)
	{
		bt->onsetDF[i] = bt->onsetDF[i+1];
	}
	
	// add new sample at the end
	bt->onsetDF[bt->onsetDFBufferSize-1] = odf_sample;
    bt->latestODF = odf_sample;
	
	// update cumulative score
	updateCumulativeScore(bt, odf_sample);
	
	// if we are halfway between beats
	if (bt->m0 == 0) {
		predictBeat(bt);
	}
	
	// if we are at a beat
	if (bt->beatCounter == 0) {
		bt->beatDueInFrame = true;	// indicate a beat should be output
		
		// recalculate the tempo
		resampleOnsetDetectionFunction(bt);
		calculateTempo(bt);
	}
}

void btrack_set_bpm(struct btrack * bt, double bpm){
	/////////// TEMPO INDICATION RESET //////////////////
	
	// firstly make sure tempo is between 80 and 160 bpm..
	while (bpm > 160) {
		bpm = bpm/2;
	}
	while (bpm < 80) {
		bpm = bpm * 2;
	}
		
	// convert tempo from bpm value to integer index of tempo probability 
	int tempo_index = (int) round((bpm - 80)/2);
	
	// now set previous tempo observations to zero
	for (int i=0;i < 41;i++) {
		bt->prevDelta[i] = 0;
	}
	
	// set desired tempo index to 1
	bt->prevDelta[tempo_index] = 1;
	
	/////////// CUMULATIVE SCORE ARTIFICAL TEMPO UPDATE //////////////////
	
	// calculate new beat period
	int new_bperiod = (int) round(60/((((double) bt->hopSize)/44100)*bpm));
	
	int bcounter = 1;
	// initialise df_buffer to zeros
	for (int i = (bt->onsetDFBufferSize-1);i >= 0;i--) {
		if (bcounter == 1) {
			bt->cumulativeScore[i] = 150;
			bt->onsetDF[i] = 150;
		} else {
			bt->cumulativeScore[i] = 10;
			bt->onsetDF[i] = 10;
		}
		
		bcounter++;
		
		if (bcounter > new_bperiod) {
			bcounter = 1;
		}
	}
	
	/////////// INDICATE THAT THIS IS A BEAT //////////////////
	
	// beat is now
	bt->beatCounter = 0;
	
	// offbeat is half of new beat period away
	bt->m0 = (int) round(((double) new_bperiod)/2);
}

void btrack_fix_bpm(struct btrack * bt, double bpm){
	// firstly make sure tempo is between 80 and 160 bpm..
	while (bpm > 160) {
		bpm = bpm/2;
	}
	while (bpm < 80) {
		bpm = bpm * 2;
	}
	
	// convert tempo from bpm value to integer index of tempo probability 
	int tempo_index = (int) round((bpm - 80)/2);
	
	// now set previous fixed previous tempo observation values to zero
	for (int i=0;i < 41;i++) {
		bt->prevDeltaFixed[i] = 0;
	}
	
	// set desired tempo index to 1
	bt->prevDeltaFixed[tempo_index] = 1;
		
	// set the tempo fix flag
	bt->tempoFixed = true;
}

void btrack_nofix_bpm(struct btrack * bt){
	// set the tempo fix flag
	bt->tempoFixed = false;
}

void btrack_set_hop_size(struct btrack * bt, int hop_size){
    bt->hopSize = hop_size;
	bt->onsetDFBufferSize = (512*512)/bt->hopSize;		// calculate df buffer size

    free(bt->onsetDF);
    bt->onsetDF = malloc(bt->onsetDFBufferSize * sizeof(double));
    BTRACK_ASSERT(bt->onsetDF);

    free(bt->cumulativeScore);
    bt->cumulativeScore = malloc(bt->onsetDFBufferSize * sizeof(double));
    BTRACK_ASSERT(bt->cumulativeScore);

    free(bt->w1);
    bt->w1 = malloc((1 + bt->onsetDFBufferSize) * sizeof(double));
    BTRACK_ASSERT(bt->w1);

    bt->onsetDFBufferSize = (512*512)/bt->hopSize;      // calculate df buffer size
    
    bt->beatPeriod = round(60/((((double) bt->hopSize)/44100)*bt->tempo));

    // initialise df_buffer to zeros
    for (int i = 0;i < bt->onsetDFBufferSize;i++) {
        bt->onsetDF[i] = 0;
        bt->cumulativeScore[i] = 0;
        if ((i %  ((int) round(bt->beatPeriod))) == 0) {
            bt->onsetDF[i] = 1;
        }
    }
}

//=======================================================================

static void resampleOnsetDetectionFunction(struct btrack * bt) {
	float output[512];
    float * input = malloc(bt->onsetDFBufferSize * sizeof(float));
    BTRACK_ASSERT(input);
    
    for (int i = 0;i < bt->onsetDFBufferSize;i++) {
        input[i] = (float) bt->onsetDF[i];
    }
		
	double src_ratio = 512.0/((double) bt->onsetDFBufferSize);
	int BUFFER_LEN = bt->onsetDFBufferSize;
	int output_len;
	SRC_DATA	src_data ;
	
	//output_len = (int) floor (((double) BUFFER_LEN) * src_ratio) ;
	output_len = 512;
	
	src_data.data_in = input;
	src_data.input_frames = BUFFER_LEN;
	
	src_data.src_ratio = src_ratio;
	
	src_data.data_out = output;
	src_data.output_frames = output_len;
	
	src_simple (&src_data, SRC_SINC_BEST_QUALITY, 1);
			
	for (int i = 0;i < output_len;i++) {
		bt->resampledOnsetDF[i] = (double) src_data.data_out[i];
	}
    free(input);
}

static void calculateTempo(struct btrack * bt){
    double aux_buffer[512];
	// adaptive threshold on input
	adaptiveThreshold(bt->resampledOnsetDF,512, aux_buffer);
		
	// calculate auto-correlation function of detection function
	calculateBalancedACF(bt->resampledOnsetDF, bt->acf);
	
	// calculate output of comb filterbank
	calculateOutputOfCombFilterBank(bt);
	
	// adaptive threshold on rcf
	adaptiveThreshold(bt->combFilterBankOutput,128, aux_buffer);

	
	int t_index;
	int t_index2;
	// calculate tempo observation vector from beat period observation vector
	for (int i = 0;i < 41;i++) {
		t_index = (int) round(bt->tempoToLagFactor / ((double) ((2*i)+80)));
		t_index2 = (int) round(bt->tempoToLagFactor / ((double) ((4*i)+160)));

		
		bt->tempoObservationVector[i] = bt->combFilterBankOutput[t_index-1] + bt->combFilterBankOutput[t_index2-1];
	}
	
	
	double maxval;
	double maxind;
	double curval;
	
	// if tempo is fixed then always use a fixed set of tempi as the previous observation probability function
	if(bt->tempoFixed) {
		for (int k = 0;k < 41;k++) {
			bt->prevDelta[k] = bt->prevDeltaFixed[k];
		}
	}
		
	for(int j=0;j < 41;j++) {
		maxval = -1;
		for (int i = 0;i < 41;i++) {
			curval = bt->prevDelta[i]*bt->tempoTransitionMatrix[i][j];
			
			if (curval > maxval) {
				maxval = curval;
			}
		}
		
		bt->delta[j] = maxval*bt->tempoObservationVector[j];
	}
	

	normaliseArray(bt->delta,41);
	
	maxind = -1;
	maxval = -1;
	
	for (int j=0;j < 41;j++) {
		if (bt->delta[j] > maxval) {
			maxval = bt->delta[j];
			maxind = j;
		}
		bt->prevDelta[j] = bt->delta[j];
	}
	
	bt->beatPeriod = round((60.0*44100.0)/(((2*maxind)+80)*((double) bt->hopSize)));
	if (bt->beatPeriod > 0) {
		bt->estimatedTempo = 60.0/((((double) bt->hopSize) / 44100.0)*bt->beatPeriod);
	}
}

static void adaptiveThreshold(double * x, int N, double * aux_buffer) {
	int i = 0;
	int k,t = 0;
	double * x_thresh = aux_buffer;
	
	int p_post = 7;
	int p_pre = 8;
	
	t = BTRACK_MIN(N,p_post);	// what is smaller, p_post of df size. This is to avoid accessing outside of arrays
	
	// find threshold for first 't' samples, where a full average cannot be computed yet 
	for (i = 0;i <= t;i++) {	
		k = BTRACK_MIN((i+p_pre),N);
		x_thresh[i] = calculateMeanOfArray(x,1,k);
	}
	// find threshold for bulk of samples across a moving average from [i-p_pre,i+p_post]
	for (i = t+1;i < N-p_post;i++) {
		x_thresh[i] = calculateMeanOfArray(x,i-p_pre,i+p_post);
	}
	// for last few samples calculate threshold, again, not enough samples to do as above
	for (i = N-p_post;i < N;i++) {
		k = BTRACK_MAX((i-p_post),1);
		x_thresh[i] = calculateMeanOfArray(x,k,N);
	}
	
	// subtract the threshold from the detection function and check that it is not less than 0
	for (i = 0;i < N;i++) {
		x[i] = x[i] - x_thresh[i];
		if (x[i] < 0) {
			x[i] = 0;
		}
	}
}

static void calculateOutputOfCombFilterBank(struct btrack * bt) {
	int numelem;
	
	for (int i = 0;i < 128;i++) {
		bt->combFilterBankOutput[i] = 0;
	}
	
	numelem = 4;
	
    // max beat period
	for (int i = 2;i <= 127;i++) {
        // number of comb elements
		for (int a = 1;a <= numelem;a++) {
            // general state using normalisation of comb elements
			for (int b = 1-a;b <= a-1;b++) {
				bt->combFilterBankOutput[i-1] = bt->combFilterBankOutput[i-1] + (bt->acf[(a*i+b)-1]*bt->weightingVector[i-1])/(2*a-1);	// calculate value for comb filter row
			}
		}
	}
}

static void calculateBalancedACF(double *onsetDetectionFunction, double * acf) {
	int l, n = 0;
	double sum, tmp;
	
	// for l lags from 0-511
	for (l = 0;l < 512;l++) {
		sum = 0;	
		
		// for n samples from 0 - (512-lag)
		for (n = 0;n < (512-l);n++) {
			tmp = onsetDetectionFunction[n] * onsetDetectionFunction[n+l];	// multiply current sample n by sample (n+l)
			sum = sum + tmp;	// add to sum
		}
		
		acf[l] = sum / (512-l);		// weight by number of mults and add to acf buffer
	}
}

static double calculateMeanOfArray(double *array, int startIndex,int endIndex) {
	int i;
	double sum = 0;

    int length = endIndex - startIndex;
	
	// find sum
	for (i = startIndex;i < endIndex;i++) {
		sum = sum + array[i];
	}
	
    if (length > 0) {
        return sum / length;	// average and return
    } else {
        return 0;
    }
}

static void normaliseArray(double * array, int N){
	double sum = 0;
	
	for (int i = 0;i < N;i++) {
		if (array[i] > 0) {
			sum = sum + array[i];
		}
	}
	
	if (sum > 0) {
		for (int i = 0;i < N;i++) {
			array[i] = array[i] / sum;
		}
	}
}

static void updateCumulativeScore(struct btrack * bt, double odfSample) {	 
	int start, end, winsize;
	double max, maxw, min;
	
	start = bt->onsetDFBufferSize - round(2*bt->beatPeriod);
	end = bt->onsetDFBufferSize - round(bt->beatPeriod/2);
	winsize = end-start+1;
	
	double v = -2*bt->beatPeriod;
	double wcumscore;
	
	
	// create window
	for (int i = 0;i < winsize;i++) {
		bt->w1[i] = exp((-1*pow(bt->tightness*log(-v/bt->beatPeriod),2))/2);
		v = v+1;
	}	
	
	// calculate new cumulative score value
	maxw = 0;
    max = 0;
    min = bt->cumulativeScore[start];
	int n = 0;
	for (int i=start;i <= end;i++) {
        wcumscore = bt->cumulativeScore[i]*bt->w1[n];
        if (wcumscore > maxw) {
            maxw = wcumscore;
        }
        if (bt->cumulativeScore[i] < min) {
            min = bt->cumulativeScore[i];
        }
        if (bt->cumulativeScore[i] > max) {
            max = bt->cumulativeScore[i];
        }
		n++;
	}
	
	
	// shift cumulative score back one
	for (int i = 0;i < (bt->onsetDFBufferSize-1);i++) {
		bt->cumulativeScore[i] = bt->cumulativeScore[i+1];
	}
	
	// add new value to cumulative score
	bt->cumulativeScore[bt->onsetDFBufferSize-1] = ((1-bt->alpha)*odfSample) + (bt->alpha*maxw);
	bt->latestCumulativeScoreValue = bt->cumulativeScore[bt->onsetDFBufferSize-1];
    bt->latestConfidence = 1.0 - (min / max);
}

static void predictBeat(struct btrack * bt){
	int windowSize = (int) bt->beatPeriod;
	double futureCumulativeScore[bt->onsetDFBufferSize + windowSize + 1];
    double w2[windowSize + 1];
    //double * futureCumulativeScore = malloc((bt->onsetDFBufferSize + windowSize + 1) * sizeof(double));
    //BTRACK_ASSERT(futureCumulativeScore);
    //double * w2 = malloc((windowSize + 1) * sizeof(double));
    //BTRACK_ASSERT(w2);

	// copy cumscore to first part of fcumscore
	for (int i = 0;i < bt->onsetDFBufferSize;i++) {
		futureCumulativeScore[i] = bt->cumulativeScore[i];
	}
	
	// create future window
	double v = 1;
	for (int i = 0;i < windowSize;i++) {
		w2[i] = exp((-1*pow((v - (bt->beatPeriod/2)),2))   /  (2*pow((bt->beatPeriod/2) ,2)));
		v++;
	}
	
	// create past window
	v = -2*bt->beatPeriod;
	int start = bt->onsetDFBufferSize - round(2*bt->beatPeriod);
	int end = bt->onsetDFBufferSize - round(bt->beatPeriod/2);
	int pastwinsize = end-start+1;

	for (int i = 0;i < pastwinsize;i++) {
		bt->w1[i] = exp((-1*pow(bt->tightness*log(-v/bt->beatPeriod),2))/2);
		v = v+1;
	}

	// calculate future cumulative score
	double max;
	int n;
	double wcumscore;
	for (int i = bt->onsetDFBufferSize;i < (bt->onsetDFBufferSize + windowSize); i++) {
		start = i - round(2*bt->beatPeriod);
		end = i - round(bt->beatPeriod/2);
		
		max = 0;
		n = 0;
		for (int k=start;k <= end;k++) {
			wcumscore = futureCumulativeScore[k]*bt->w1[n];
			
			if (wcumscore > max) {
				max = wcumscore;
			}
			n++;
		}
		
		futureCumulativeScore[i] = max;
	}
	
	
	// predict beat
	max = 0;
	n = 0;
	
	for (int i = bt->onsetDFBufferSize;i < (bt->onsetDFBufferSize+windowSize);i++) {
		wcumscore = futureCumulativeScore[i]*w2[n];
		if (wcumscore > max) {
			max = wcumscore;
			bt->beatCounter = n;
		}	
		n++;
	}
		
	// set next prediction time
	bt->m0 = bt->beatCounter+round(bt->beatPeriod/2);

    //free(futureCumulativeScore);
    //free(w2);
}
