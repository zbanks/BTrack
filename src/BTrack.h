//=======================================================================
/** @file BTrack.h
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

#ifndef __BTRACK_H
#define __BTRACK_H

#include "OnsetDetectionFunction.h"

//=======================================================================
/** The main beat tracking class and the interface to the BTrack
 * beat tracking algorithm. The algorithm can process either
 * audio frames or onset detection function samples and also
 * contains some static functions for calculating beat times in seconds
 */
struct btrack {
    struct odf odf;
    int frameSize;
    double invSampleRate;
    double * onsetDF;                       /**< to hold onset detection function */
    double * cumulativeScore;               /**< to hold cumulative score */
    double * w1;

    double resampledOnsetDF[512];           /**< to hold resampled detection function */
    double acf[512];                        /**<  to hold autocorrelation function */
    double weightingVector[128];            /**<  to hold weighting vector */
    double combFilterBankOutput[128];       /**<  to hold comb filter output */
    double tempoObservationVector[41];      /**<  to hold tempo version of comb filter output */
	
    double delta[41];                       /**<  to hold final tempo candidate array */
    double prevDelta[41];                   /**<  previous delta */
    double prevDeltaFixed[41];              /**<  fixed tempo version of previous delta */
	
    double tempoTransitionMatrix[41][41];   /**<  tempo transition matrix */
    
	//=======================================================================
    // parameters
    
    double tightness;                       /**< the tightness of the weighting used to calculate cumulative score */
    double alpha;                           /**< the mix between the current detection function sample and the cumulative score's "momentum" */
    double beatPeriod;                      /**< the beat period, in detection function samples */
    double tempo;                           /**< the tempo in beats per minute */
    double estimatedTempo;                  /**< the current tempo estimation being used by the algorithm */
    double latestCumulativeScoreValue;      /**< holds the latest value of the cumulative score function */
    double latestODF;                       /**< holds the latest value of the onset detection function*/
    double latestConfidence;                /**< holds the latest confidence value, the ratio between max score & min score in the last beat */
    double tempoToLagFactor;                /**< factor for converting between lag and tempo */
    int m0;                                 /**< indicates when the next point to predict the next beat is */
    int beatCounter;                        /**< keeps track of when the next beat is - will be zero when the beat is due, and is set elsewhere in the algorithm to be positive once a beat prediction is made */
    int hopSize;                            /**< the hop size being used by the algorithm */
    int onsetDFBufferSize;                  /**< the onset detection function buffer size */
    int tempoFixed;                         /**< indicates whether the tempo should be fixed or not */
    int beatDueInFrame;                     /**< indicates whether a beat is due in the current frame */

};

/** Constructor taking both hopSize and frameSize
* @param hop_size the hop size in audio samples
*   - hop_size determines how often the state is updated
* @param frame_size the frame size in audio samples
*   - frame_size samples are used to calculate the ODF
* Constraint: hop_size <= frame_size
*/
int btrack_init(struct btrack * bt, int hop_size, int frame_size, int sample_rate);
void btrack_del(struct btrack * bt);

void btrack_process_audio_frame(struct btrack * bt, const btrack_chunk_t * frame);
void btrack_process_fft_frame(struct btrack * bt, const btrack_chunk_t * fft_frame);
void btrack_process_odf_sample(struct btrack * bt, double odf_sample);

int btrack_beat_due_in_current_frame(const struct btrack * bt);
double btrack_get_bpm(const struct btrack * bt);
double btrack_get_latest_score(const struct btrack * bt);
double btrack_get_latest_odf(const struct btrack * bt);
double btrack_get_latest_confidence(const struct btrack * bt);
int btrack_get_frames_until_beat(const struct btrack * bt);
double btrack_get_time_until_beat(const struct btrack * bt);

// Reset the internal probability state to a given bpm
// Useful for giving btrack a "hint" about the current bpm
// Over time, the tempo can drift
void btrack_set_bpm(struct btrack * bt, double bpm);

// Lock the tempo to the given bpm
// The tempo may have small deviations in order to correct for phase
void btrack_fix_bpm(struct btrack * bt, double bpm);
void btrack_nofix_bpm(struct btrack * bt);

void btrack_set_hop_size(struct btrack * bt, int hop_size);

#endif
