//=======================================================================
/** @file OnsetDetectionFunction.h
 *  @brief A class for calculating onset detection functions
 *  @author Adam Stark
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

#ifndef __ONSETDETECTIONFUNCTION_H
#define __ONSETDETECTIONFUNCTION_H

#include "fftw3.h"

//=======================================================================
/** The type of onset detection function to calculate */
enum OnsetDetectionFunctionType
{
    EnergyEnvelope,
    EnergyDifference,
    SpectralDifference,
    SpectralDifferenceHWR,
    PhaseDeviation,
    ComplexSpectralDifference,
    ComplexSpectralDifferenceHWR,
    HighFrequencyContent,
    HighFrequencySpectralDifference,
    HighFrequencySpectralDifferenceHWR
};

//=======================================================================
/** The type of window to use when calculating onset detection function samples */
enum WindowType
{
    RectangularWindow,
    HanningWindow,
    HammingWindow,
    BlackmanWindow,
    TukeyWindow
};

//=======================================================================
/** A class for calculating onset detection functions. */
struct odf {
	int frameSize;						/**< audio framesize */
	int hopSize;						/**< audio hopsize */
	enum OnsetDetectionFunctionType onsetDetectionFunctionType;		/**< type of detection function */
    enum WindowType windowType;                     /**< type of window used in calculations */
	
	fftw_plan p;						/**< fftw plan */
	fftw_complex *complexIn;			/**< to hold complex fft values for input */
	fftw_complex *complexOut;			/**< to hold complex fft values for output */
	
	int initialised;					/**< flag indicating whether buffers and FFT plans are initialised */

    double * frame;                     /**< audio frame */
    double * window;                    /**< window */
	
	double prevEnergySum;				/**< to hold the previous energy sum value */
	
    double * magSpec;                   /**< magnitude spectrum */
    double * prevMagSpec;               /**< previous magnitude spectrum */
	
    double * phase;                     /**< FFT phase values */
    double * prevPhase;                 /**< previous phase values */
    double * revPhase2;                 /**< second order previous phase values */

};

/** Constructor 
* @param hopSize_ the hop size in audio samples
* @param frameSize_ the frame size in audio samples
* @param onsetDetectionFunctionType_ the type of onset detection function to use - (see OnsetDetectionFunctionType)
* @param windowType the type of window to use (see WindowType)
*/
int odf_init(struct odf * odf, int hop_size, int frame_size, enum OnsetDetectionFunctionType odf_type, enum WindowType window_type);
void odf_del(struct odf * odf);

double odf_calculate_sample(struct odf * odf, double * buffer);
void odf_set_type(struct odf * odf, enum OnsetDetectionFunctionType type);

#endif
