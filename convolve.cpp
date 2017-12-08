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
    float*         dataArray;
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

    input.close();

    return waveFile;
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
    cout << "--------------------------" << endl;
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

//convolve skeleton function gotten from professor Manzara convolution demo program
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
    printWAVEdetails(inputWave);
    Wave IRWave = readWAVE(IRfile);
    printWAVEdetails(IRWave);

    return 0;
}
