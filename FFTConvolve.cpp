#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>

using namespace std;

using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;

//Code for wav file struct gotten from user kory @ https://stackoverflow.com/questions/13660777/c-reading-the-data-part-of-a-wav-file
typedef struct  WAV
{
    /* RIFF Chunk Descriptor */
    uint8_t         RIFF[4];        // RIFF Header Magic header
    uint32_t        ChunkSize;      // RIFF Chunk Size
    uint8_t         WAVE[4];        // WAVE Header
    /* "fmt" sub-chunk */
    uint8_t         fmt[4];         // FMT header
    uint32_t        Subchunk1Size;  // Size of the fmt chunk
    uint16_t        AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
    uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Sterio
    uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
    uint32_t        bytesPerSec;    // bytes per second
    uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    uint16_t        bitsPerSample;  // Number of bits per sample

    uint8_t*         extraData;
    /* "data" sub-chunk */
    uint8_t         Subchunk2ID[4]; // "data"  string
    uint32_t        Subchunk2Size;  // Sampled data length

    //actual data
    float*          dataArray;
    int             dataLength;
} Wave;

// Function prototypes
int getFileSize(FILE* inFile);

Wave readWAVE(const char* filePath){

    Wave waveFile;
    ifstream input;

    input.open(filePath, ios::binary);

    input.read((char*)&waveFile.RIFF, 4);
    input.read((char*)&waveFile.ChunkSize, 4);
    input.read((char*)&waveFile.WAVE, 4);
    input.read((char*)&waveFile.fmt, 4);
    input.read((char*)&waveFile.Subchunk1Size, 4);
    input.read((char*)&waveFile.AudioFormat, 2);
    input.read((char*)&waveFile.NumOfChan, 2);
    input.read((char*)&waveFile.SamplesPerSec, 4);
    input.read((char*)&waveFile.bytesPerSec, 4);
    input.read((char*)&waveFile.blockAlign, 2);
    input.read((char*)&waveFile.bitsPerSample, 2);
    waveFile.extraData = new uint8_t[waveFile.Subchunk1Size - 16];
    input.read((char*)&waveFile.extraData, waveFile.Subchunk1Size-16);
    input.read((char*)&waveFile.Subchunk2ID, 4);
    input.read((char*)&waveFile.Subchunk2Size, 4);
    waveFile.dataArray = new float[waveFile.Subchunk2Size / (waveFile.bitsPerSample / 8)];
    int16_t sample;
    for(int i = 0; i < waveFile.Subchunk2Size / (waveFile.bitsPerSample / 8); i++){
        input.read((char*)&sample, 2);
        float converted = (float) sample / (float)INT16_MAX;
        if(converted < -1.0){
            converted = -1.0;
        }
        waveFile.dataArray[i] = converted;
    }
    waveFile.dataLength = (waveFile.Subchunk2Size / (waveFile.bitsPerSample/8));

    //strip one channel if stereo
    if (waveFile.NumOfChan == 2){
        float* strip = new float[waveFile.dataLength/2];
        for (int i = 0; i < waveFile.dataLength/2; i++) {
            strip[i] = waveFile.dataArray[i*2];
        }
        waveFile.dataArray = strip;
        waveFile.dataLength = waveFile.dataLength/2;
    }

    input.close();

    return waveFile;
}

void writeWAVE(string outputFile, Wave original, float data[], int dataLength) {

    char* riff = new char[4]{(char)original.RIFF[0], (char)original.RIFF[1], (char)original.RIFF[2], (char)original.RIFF[3]};       //"RIFF"
    int Subchunk2Size = dataLength * 2; //Size of sound sample data in bytes
    int chunkSize = Subchunk2Size + 36; //Size of remaining file in bytes (36 + Subchunk2Size)
    char* format = new char[4]{(char)original.WAVE[0], (char)original.WAVE[1], (char)original.WAVE[2], (char)original.WAVE[3]};     //"WAVE"
    char* subChunk1ID = new char[4]{(char)original.fmt[0], (char)original.fmt[1], (char)original.fmt[2], (char)original.fmt[3]};    //"fmt" space after t is needed
    int subChunk1Size = 16; // 16
    short audioFmt = 1; //1 = PCM
    short numChannels = 1;  //1 = mono 2 = stereo
    int byteRate = (int)original.SamplesPerSec * numChannels * original.bitsPerSample/8;    //SamepleRate * NumChannels * BitsPerSample/8
    int blockAlign = numChannels * (original.bitsPerSample / 8);    //NumChannels * BitsPerSample/8
    char* subChunk2ID = new char[4]{(char)original.Subchunk2ID[0], (char)original.Subchunk2ID[1], (char)original.Subchunk2ID[2], (char)original.Subchunk2ID[3]};    //"data"

    //write
    ofstream output;
    output.open(outputFile, ios::binary | ios::out);
    output.write(riff, 4);                  // RIFF
    output.write((char*)&chunkSize, 4);     //ChunkSize
    output.write(format, 4);                //Format
    output.write(subChunk1ID, 4);           //subChunk1ID
    output.write((char*)&subChunk1Size, 4); //subChunk1Size
    output.write((char*)&audioFmt, 2);      //AudioFormat
    output.write((char*)&numChannels, 2);   //NumChannels
    output.write((char*)&original.SamplesPerSec, 4);    //SampleRate
    output.write((char*)&byteRate, 4);      //ByteRate
    output.write((char*)&blockAlign, 2);    //BlockAlign
    output.write((char*)&original.bitsPerSample, 2);    //BitsPerSample
    output.write(subChunk2ID, 4);           //Subchunk2ID
    output.write((char*)&Subchunk2Size, 4); //Subchunk2Size
    //Data
    int16_t dataToWrite;
    for(int i = 0; i < dataLength; i++){
        dataToWrite = (int16_t) (data[i] * INT16_MAX);
        output.write((char*)&dataToWrite, 2);
    }
    output.close();
}

//Code for wav file printout format gotten from user kory @ https://stackoverflow.com/questions/13660777/c-reading-the-data-part-of-a-wav-file
void printWAVEdetails(Wave wave){

    cout << "RIFF header                :" << wave.RIFF[0] << wave.RIFF[1] << wave.RIFF[2] << wave.RIFF[3] << endl;
    cout << "WAVE header                :" << wave.WAVE[0] << wave.WAVE[1] << wave.WAVE[2] << wave.WAVE[3] << endl;
    cout << "FMT                        :" << wave.fmt[0] << wave.fmt[1] << wave.fmt[2] << wave.fmt[3] << endl;
    cout << "Data size                  :" << wave.ChunkSize << endl;

    // Display the sampling Rate from the header
    cout << "Sampling Rate              :" << wave.SamplesPerSec << endl;
    cout << "Number of bits used        :" << wave.bitsPerSample << endl;
    cout << "Number of channels         :" << wave.NumOfChan << endl;
    cout << "Number of bytes per second :" << wave.bytesPerSec << endl;
    cout << "Data length                :" << wave.Subchunk2Size << endl;
    cout << "Audio Format               :" << wave.AudioFormat << endl;
    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM

    cout << "Block align                :" << wave.blockAlign << endl;
    cout << "Data string                :" << wave.Subchunk2ID[0] << wave.Subchunk2ID[1] << wave.Subchunk2ID[2] << wave.Subchunk2ID[3] << endl;
    cout << "----------------------------------" << endl;
}

// find the file size
int getFileSize(FILE* inFile)
{
    int fileSize = 0;
    fseek(inFile, 0, SEEK_END);

    fileSize = ftell(inFile);

    fseek(inFile, 0, SEEK_SET);
    return fileSize;
}

//  From 501 class handout Dr Manzara
//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

//convolve function gotten from professor Manzara convolution demo program
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
  int n, m;

  /*  Make sure the output buffer is the right size: P = N + M - 1  */
  if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    return;
  }

  /*  Clear the output buffer y[] to all zero values  */
  for (n = 0; n < P; n++)
    y[n] = 0.0;

  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */
  for (n = 0; n < N; n++) {
    /*  Inner loop:  process x[n] with each sample of h[]  */
    for (m = 0; m < M; m++)
      y[n+m] += x[n] * h[m];
  }
}

int main(int argc, char* argv[])
{
    const char* inputfile;
    const char* IRfile;
    const char* outputfile;

    if (argc < 4) {
        cout << "Usage: " << argv[0] << " inputfile IRfile outputfile\n";
        exit(0);
    }
    else {
        inputfile = argv[1];
        cout << "input name: " << inputfile << endl;
        IRfile = argv[2];
        cout << "IRfile name: " << inputfile << endl;
        outputfile = argv[3];
        cout << "output name: " << inputfile << endl << endl;
    }
    Wave inputWave = readWAVE(inputfile);
    cout << "inputfile details\n---" << endl;
    printWAVEdetails(inputWave);
    Wave IRWave = readWAVE(IRfile);
    cout << "IRfile details\n---" << endl;
    printWAVEdetails(IRWave);

    int convolvedDataLength = inputWave.dataLength + IRWave.dataLength - 1;
    float* convolvedData = new float[convolvedDataLength];

    cout << "Convolving..." << endl;
    convolve(inputWave.dataArray, inputWave.dataLength, IRWave.dataArray, IRWave.dataLength, convolvedData, convolvedDataLength);
    writeWAVE(outputfile, inputWave, convolvedData, convolvedDataLength);
    // writeWAVE(outputfile, inputWave, inputWave.dataArray, inputWave.dataLength); //test writing without convolving
    cout << "Done!" << endl;

    Wave outputWave = readWAVE(outputfile);
    cout << "outputfile details\n---" << endl;
    printWAVEdetails(outputWave);

    cout << "Convolution reverb has been completed on " << inputfile << " with impulse response " << IRfile << " and file " << outputfile << " has been created" << endl;

    return 0;
}
