package main

import (
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"os"

	"github.com/hajimehoshi/ebiten/audio"
	"github.com/hajimehoshi/ebiten/v2/audio/wav"
	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func lowPassFilterLFEContinuous(lfe []complex128, sampleRate float64, cutoffFreq float64) {
	M := len(lfe)    // Nombre de coefficients fréquentiels
	N := 2 * (M - 1) // Longueur du signal temporel
	freqResolution := sampleRate / float64(N)

	// Paramètre de pente pour la décroissance (plus tau est grand, plus la transition est douce)
	tau := 0.7 * cutoffFreq

	fmt.Printf("freqResolution=%v Number_M_of_freq_coeff=%v\n", freqResolution, M)

	// Appliquer une fonction continue dans le domaine fréquentiel
	for i := 0; i < M; i++ {
		// Fréquence correspondante à l'indice i
		freq := float64(i) * freqResolution

		// Fonction continue : atténuation exponentielle après cutoffFreq
		if freq > cutoffFreq {
			attenuation := math.Exp(-((freq - cutoffFreq) / tau))
			lfe[i] = lfe[i] * complex(attenuation, 0)
		}
		// Les fréquences <= cutoffFreq restent inchangées (gain de 1)
	}
}

func lowPassFilterLFE(lfe []complex128, sampleRate float64, cutoffFreq float64) {
	M := len(lfe)    // Nombre de coefficients fréquentiels
	N := 2 * (M - 1) // Longueur du signal temporel
	freqResolution := sampleRate / float64(N)
	cutoffIndex := int(cutoffFreq/freqResolution) + 1

	fmt.Printf("freqResolution=%v cutoffIndex=%v Number_M_of_freq_coeff=%v\n", freqResolution, cutoffIndex, M)

	for i := cutoffIndex; i < M; i++ {
		lfe[i] = complex(0, 0) // Supprime les hautes fréquences positives
	}
}

func maxAbs(data []float64) float64 {
	max := 0.0
	for _, v := range data {
		abs := math.Abs(v)
		if abs > max {
			max = abs
		}
	}
	return max
}

func normalize(left *[]float64, right *[]float64) {
	maxVal := math.Max(maxAbs(*left), maxAbs(*right))
	if maxVal > 1 {
		for i := range *left {
			(*left)[i] /= maxVal
			(*right)[i] /= maxVal
		}
	}
}

func normalizeSingle(channel *[]float64) {
	maxVal := maxAbs(*channel)
	if maxVal > 1 {
		for i := range *channel {
			(*channel)[i] /= maxVal
		}
	}
}

func plotSignal(x []float64, xFiltered []float64, filename string) {
	p := plot.New()
	p.Title.Text = "Before and After Filter Signal"
	p.X.Label.Text = "Time"
	p.Y.Label.Text = "Amplitude"

	ptsOriginal := make(plotter.XYs, len(x))
	ptsFiltered := make(plotter.XYs, len(xFiltered))

	for i := range x {
		ptsOriginal[i].X = float64(i)
		ptsOriginal[i].Y = x[i]
		ptsFiltered[i].X = float64(i)
		ptsFiltered[i].Y = xFiltered[i]
	}

	// Signal original en vert (couleur par défaut)
	lineOriginal, err := plotter.NewLine(ptsOriginal)
	if err != nil {
		panic(err)
	}
	lineOriginal.Color = plotutil.Color(0) // Bleu

	// Signal filtré en rouge
	lineFiltered, err := plotter.NewLine(ptsFiltered)
	if err != nil {
		panic(err)
	}
	lineFiltered.Color = plotutil.Color(1) // Rouge
	lineFiltered.LineStyle.Width = vg.Points(2)

	p.Add(lineOriginal, lineFiltered)
	p.Legend.Add("Original", lineOriginal)
	p.Legend.Add("Filtered", lineFiltered)

	if err := p.Save(8*vg.Inch, 4*vg.Inch, filename); err != nil {
		panic(err)
	}
}

func plotSpectrum(X []complex128, XFiltered []complex128, sampleRate float64, filename string) {
	p := plot.New()
	p.Title.Text = "Signal Spectrum Before and After Filtering"
	p.X.Label.Text = "Frequence (Hz)"
	p.Y.Label.Text = "Amplitude"

	N := len(X) // Nombre de fréquences analysées
	freqResolution := sampleRate / float64(2*(N-1))

	ptsOriginal := make(plotter.XYs, N)
	ptsFiltered := make(plotter.XYs, N)

	for k := 0; k < N; k++ {
		freq := float64(k) * freqResolution
		ptsOriginal[k].X = freq
		ptsOriginal[k].Y = cmplx.Abs(X[k]) // Module des coefficients FFT
		ptsFiltered[k].X = freq
		ptsFiltered[k].Y = cmplx.Abs(XFiltered[k])
	}

	// Spectre original en bleu
	lineOriginal, _ := plotter.NewLine(ptsOriginal)
	lineOriginal.Color = plotutil.Color(0) // Bleu

	// Spectre filtré en rouge
	lineFiltered, _ := plotter.NewLine(ptsFiltered)
	lineFiltered.Color = plotutil.Color(1) // Rouge
	lineFiltered.LineStyle.Width = vg.Points(2)

	p.Add(lineOriginal, lineFiltered)
	p.Legend.Add("Original Spectrum", lineOriginal)
	p.Legend.Add("Filtered Spectrum", lineFiltered)

	_ = p.Save(8*vg.Inch, 4*vg.Inch, filename)
}

func createWAVHeader(sampleRate int) []byte {

	// Write WAV header manually
	// >> is used to perform right bit shift
	// sampleRate >> 8 shifts the bits 8 positions to the right, effectively dividing by 256 (2^8) and getting the second byte of the number.
	// sampleRate >> 16 shifts 16 positions, giving the third byte.
	// sampleRate >> 24 for the fourth byte.
	// The & 0xFF operation masks out all but the least significant byte after the shift, ensuring only one byte is written.

	return []byte{
		'R', 'I', 'F', 'F', 0, 0, 0, 0, // RIFF
		'W', 'A', 'V', 'E', // WAVE
		'f', 'm', 't', ' ', 16, 0, 0, 0, // fmt
		1, 0, // Compression code (1 = PCM)
		2, 0, // Number of channels
		byte(sampleRate & 0xFF), byte((sampleRate >> 8) & 0xFF), byte((sampleRate >> 16) & 0xFF), byte((sampleRate >> 24) & 0xFF), // Sample rate
		byte((sampleRate * 4) & 0xFF), byte(((sampleRate * 4) >> 8) & 0xFF), byte(((sampleRate * 4) >> 16) & 0xFF), byte(((sampleRate * 4) >> 24) & 0xFF), // Byte rate (sampleRate * channels * bitsPerSample / 8)
		4, 0, // Block align
		16, 0, // Bits per sample
		'd', 'a', 't', 'a', 0, 0, 0, 0, // data
	}
}

func writeWaveFile(s string, sampleRate int, left []float64, right []float64) error {
	// Write back channels to new WAV file
	outFile, err := os.Create(s)
	if err != nil {
		return fmt.Errorf("error creating  WAV file: %w", err)

	}
	defer outFile.Close()

	header := createWAVHeader(sampleRate)

	// Write the header to the file
	outFile.Write(header)

	// Write audio data
	// Each audio signal  is converted from float to a 16-bit integer, then split into two bytes
	// & 0xFF gets the least significant byte.
	// >> 8 shifts to get the next byte for the 16-bit integer representation.

	for i := 0; i < len(left); i++ {
		outFile.Write([]byte{
			byte(int16(left[i]*float64(math.MaxInt16)) & 0xFF),
			byte(int16(left[i]*float64(math.MaxInt16)) >> 8),
			byte(int16(right[i]*float64(math.MaxInt16)) & 0xFF),
			byte(int16(right[i]*float64(math.MaxInt16)) >> 8),
		})
	}

	// Update header sizes
	outSize := int64(len(left) * 4) // Each sample is 4 bytes (2 channels * 2 bytes per sample)
	chunkSize := 36 + outSize       // 36 = size of the header up to data chunk
	header[4] = byte(chunkSize & 0xFF)
	header[5] = byte((chunkSize >> 8) & 0xFF)
	header[6] = byte((chunkSize >> 16) & 0xFF)
	header[7] = byte((chunkSize >> 24) & 0xFF)
	header[40] = byte(outSize & 0xFF)
	header[41] = byte((outSize >> 8) & 0xFF)
	header[42] = byte((outSize >> 16) & 0xFF)
	header[43] = byte((outSize >> 24) & 0xFF)

	// Write the updated header back to the file
	outFile.Seek(0, 0)
	outFile.Write(header)

	return nil
}

func readWaveFile(s string) ([]float64, []float64, int, error) {
	f, err := os.Open(s)
	if err != nil {
		return nil, nil, 0, fmt.Errorf("error opening WAV file: %w", err)
	}
	defer f.Close()

	context, err := audio.NewContext(44100) // Assuming 44.1 kHz; adjust if necessary
	if err != nil {
		return nil, nil, 0, fmt.Errorf("failed to create audio context: %w", err)
	}

	d, err := wav.DecodeWithSampleRate(context.SampleRate(), f)
	if err != nil {
		return nil, nil, 0, fmt.Errorf("error decoding WAV: %w", err)
	}

	var sampleRate int = context.SampleRate()
	err = nil

	data, err := io.ReadAll(d)
	if err != nil {
		return nil, nil, 0, fmt.Errorf("error reading WAV data: %w", err)
	}

	// Convert byte data to float64 for processing
	LT := make([]float64, len(data)/4)
	RT := make([]float64, len(data)/4)

	for i := 0; i < len(data)/4; i++ {
		LT[i] = float64(int16(data[i*4+0])|int16(data[i*4+1])<<8) / float64(math.MaxInt16)
		RT[i] = float64(int16(data[i*4+2])|int16(data[i*4+3])<<8) / float64(math.MaxInt16)
	}

	return LT, RT, sampleRate, nil
}

func FFTaudio() {
	LT, RT, sampleRate, err := readWaveFile("sqdemo1.wav")
	if err != nil {
		fmt.Printf("Failed to read: %v", err)
		return
	}

	normalize(&LT, &RT)

	N := len(LT)
	fft := fourier.NewFFT(N)
	freqLT := fft.Coefficients(nil, LT)
	freqRT := fft.Coefficients(nil, RT)
	M := len(freqLT)

	freqLTo := make([]complex128, M)
	copy(freqLTo, freqLT)
	freqRTo := make([]complex128, M)
	copy(freqRTo, freqRT)

	cutoff := 150.0

	// multiplying the  spectrum by a rectangular window function H(f) in the frequency domain.
	// lowPassFilterLFE(freqLT, float64(sampleRate), cutoff)
	// lowPassFilterLFE(freqRT, float64(sampleRate), cutoff)

	// multiplying the spectrum by a continuous function H(f) in the frequency domain.
	lowPassFilterLFEContinuous(freqLT, float64(sampleRate), cutoff)
	lowPassFilterLFEContinuous(freqRT, float64(sampleRate), cutoff)

	frontLeftTime := fft.Sequence(nil, freqLT)
	frontRightTime := fft.Sequence(nil, freqRT)

	normalize(&frontLeftTime, &frontRightTime)

	err = writeWaveFile("sqdemo1filtered.wav", sampleRate, frontLeftTime, frontRightTime)
	if err != nil {
		fmt.Printf("Failed to write output back channels: %v", err)
		return
	}

	plotSignal(LT, frontLeftTime, "filtered_audio_signal.png")
	plotSpectrum(freqLTo, freqLT, float64(sampleRate), "spectrum_audio.png")

}

func FFTSinusoid() {
	// N is the length of the input time-domain signal
	// sampleRate is the number of samples per second (in Hz) used to capture the signal. It determines the range of frequencies we can analyze.
	// The duration T of the signal in the time domain: T = N / sampleRate
	// freqResolution  = Δf = sampleRate / N
	// Shannon-Nyquist Theorem  fmax = sampleRate/2
	// Sampling Period Δt = 1 / sampleRate

	sampleRate := 750.0 //  Sampling Frequency
	N := 32             // Number of samples (power of 2 preferred)
	x := make([]float64, N)

	fmt.Print("Sinusoïd\n")

	for n := 0; n < N; n++ {
		// Signal réel
		x[n] = math.Sin(2 * math.Pi * 50 * float64(n) / sampleRate)
		x[n] = x[n] + math.Sin(math.Pi*50*float64(n)/sampleRate)
		x[n] = x[n] - math.Cos(2*math.Pi*75*float64(n)/sampleRate)
	}

	normalizeSingle(&x)

	fmt.Printf(" %v \n", x)

	// FFT optimisée pour les signaux réels
	fft := fourier.NewFFT(N)
	X := fft.Coefficients(nil, x) // La sortie a une longueur de N/2 + 1

	// Affichage des fréquences et des coefficients FFT
	for k := 0; k < len(X); k++ {
		freq := float64(k) * sampleRate / float64(N)
		fmt.Printf("Freq: %.1f Hz, FFT: %v\n", freq, X[k])
	}

	Xoriginal := make([]complex128, len(X))
	copy(Xoriginal, X)

	cutoff := 67.0
	lowPassFilterLFEContinuous(X, sampleRate, cutoff)

	for k := 0; k < len(X); k++ {
		freq := float64(k) * sampleRate / float64(N)
		fmt.Printf("Freq Filtered: %.1f Hz, cutoff %.1f Hz,FFT: %v\n", freq, cutoff, X[k])
	}

	xFiltered := fft.Sequence(nil, X)
	normalizeSingle(&xFiltered)

	fmt.Println("Signal original :\n", x)
	fmt.Println("Signal filtré :\n", xFiltered)

	plotSignal(x, xFiltered, "filtered__sinusoid_signal.png")
	plotSpectrum(Xoriginal, X, sampleRate, "spectrum_sinusoid.png")

}

func main() {
	FFTSinusoid()
	FFTaudio()
}
