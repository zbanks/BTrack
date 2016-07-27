#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sndfile.h>
#include <BTrack.h>

int main(int argc, char ** argv) {
    if (argc < 2) {
        printf("Usage: ./test file.wav\n");
        return 1;
    }

    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof sfinfo);
    SNDFILE * sndfile = sf_open(argv[1], SFM_READ, &sfinfo);
    if (sndfile == NULL) {
        printf("Unable to open file %s; %s\n", argv[1], sf_strerror(sndfile));
        return 1;
    }

    int block_size = 512;
    float block[512];

    struct btrack btrack;
    int rc = btrack_init(&btrack, block_size, block_size * 2, sfinfo.samplerate);
    if (rc < 0) return 1;

    if (sfinfo.channels != 1) {
        printf("Expected 1 channel, got %d\n", sfinfo.channels);
        return 1;
    }

    double time = 0;
    time -= 1.5 * (double) block_size / (double) sfinfo.samplerate;
    while (1) {
        int count = sf_readf_float(sndfile, block, block_size);
        if (count < 0) return 1;
        if (count != block_size) break;

        btrack_process_audio_frame(&btrack, block);
        if (btrack_beat_due_in_current_frame(&btrack))
            printf("%0.8f\n", time);
        time += (double) block_size / (double) sfinfo.samplerate;
        
    }
    //printf("Total time: %0.3f\n", time);

    btrack_del(&btrack);
    sf_close(sndfile);
    return 0;
}
