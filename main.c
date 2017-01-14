//
//  main.c
//  spectrogram
//
//  Created by Jeremy on 12/17/15.
//  Copyright (c) 2015 Jeremy. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define sRate 16000 //sample Rate
#define PI 3.14159
#define ANALYSIS_WINDOW_SIZE_IN_SMP 512
#define CEPSTRUM_ORDER 20

float *ham;

void HammingWindowCon(size_t len)
{

    size_t i;
    float radian ;
    ham = ( float *) calloc(sizeof(float), len);
    for (i=0; i<len;i++) {
        radian = 2.0 *PI*((float)i)/((float)(len-1));
        ham[i]=0.54-0.46*cos(radian);
    }
}
void HammingClean(void){
    if (ham) {
        free(ham);
    }
}

void dft(int N,float *time_real_part,float *time_imaginary_part,float *frequency_real_part,float *frequency_imaginary_part){

    int n,k;
    float tmp;
    float radian;
    for (k=0;k<N ; k++) {
        frequency_real_part[k] = 0;
        frequency_imaginary_part[k] = 0;
        for (n=0; n<N; n++) {
            radian = 2.0*PI/((float)N)*((float)k)*((float)n);
            tmp = time_real_part[n]*cos(radian)+(time_imaginary_part[n]) * sin(radian);
            frequency_real_part[k] += tmp;
            tmp = time_imaginary_part[n]* cos(radian) - (time_real_part[n])*sin(radian);
            frequency_imaginary_part[k] += tmp;
        }
        frequency_real_part[k] /= ((float)N);//scale
        frequency_imaginary_part[k] /= ((float)N);//scale
    
    }




}



int main(void) {
    // insert code here...
    char pcmfn[256] = ("jeremy.wav");//insert wav
    size_t total_len_in_byte = 0;
    size_t pcm_len_in_byte = 0;
    size_t pcm_len_in_smp = 0;
    FILE *fp = fopen(pcmfn, "rb");
    FILE *fp_save = NULL;
    FILE *fp_save_C = fopen("cepstrum.txt", "w");//create a cepstrum.txt
    FILE *fp_save_Reespec = fopen("reconstructed_spec.txt", "w");// create a reconstructed_spec.txt
    
    short *pcmbuffer = NULL;
    float *xm;
    float *magnitude;
    
    size_t n = 0;
    size_t m;
    float *time_real_part = NULL, *time_imaginary_part = NULL, *frequency_real_part = NULL, *frequency_imaginary_part = NULL, *cepstrum = NULL, *Realspectrum = NULL;
    
    size_t window_size_in_ms = 30, frame_size_in_ms = 10;/////millisecond
    size_t window_size_in_smp = sRate * window_size_in_ms / 1000;
    size_t frame_size_in_smp = sRate * frame_size_in_ms /1000;
    
    fseek(fp, 0, SEEK_END);
    total_len_in_byte = ftell(fp);
    pcm_len_in_byte = total_len_in_byte -44;
    pcm_len_in_smp = pcm_len_in_byte / 2;//short is 16 bits(2 bytes)
    fseek(fp, 44, SEEK_SET);
    
    pcmbuffer = (short *) calloc(sizeof(short), pcm_len_in_smp);
    time_real_part = (float *) calloc(sizeof(float), pcm_len_in_smp);
    time_imaginary_part = (float *) calloc(sizeof(float), pcm_len_in_smp);
    frequency_real_part = (float *) calloc(sizeof(float), pcm_len_in_smp);
    frequency_imaginary_part = (float *) calloc(sizeof(float), pcm_len_in_smp);
    cepstrum = (float *) calloc(sizeof(float), ANALYSIS_WINDOW_SIZE_IN_SMP);
    Realspectrum = (float *) calloc(sizeof(float), ANALYSIS_WINDOW_SIZE_IN_SMP);
    xm = (float *) calloc(sizeof(float), ANALYSIS_WINDOW_SIZE_IN_SMP);
    magnitude = (float *) calloc(sizeof(float), ANALYSIS_WINDOW_SIZE_IN_SMP);
    
    HammingWindowCon(window_size_in_smp);
    fread(pcmbuffer, sizeof(short), pcm_len_in_smp, fp);//read
    fprintf(stdout, "time.domain\n");
    for(n = 0; n<pcm_len_in_smp;n++){
        time_real_part[n] = (float)(pcmbuffer[n]);
        fprintf(stdout, "%d\n",pcmbuffer[n]);

    }
    
    fclose(fp);

    n= 0;
    fp_save = fopen("spectrogram.txt", "w");
    while ( n < pcm_len_in_smp) {
        memset(xm, 0, sizeof(float)* ANALYSIS_WINDOW_SIZE_IN_SMP);
        memcpy(xm, time_real_part + n, sizeof(float)* window_size_in_smp);
        
        for ( m =0; m < window_size_in_smp; m++) {
            xm[m] = xm[m] *ham[m];
        }
        dft(ANALYSIS_WINDOW_SIZE_IN_SMP, xm, time_imaginary_part, frequency_real_part, frequency_imaginary_part);
        for (m=0; m< (ANALYSIS_WINDOW_SIZE_IN_SMP); m++) {
            magnitude[m] = 10.0 * log10(frequency_real_part[m] * frequency_real_part[m] + frequency_imaginary_part[m] * frequency_imaginary_part[m]);
            fprintf(fp_save,"%.15f\t",magnitude[m]);
            fflush(fp_save);
        }
        
        fprintf(fp_save, "\n");
        dft(ANALYSIS_WINDOW_SIZE_IN_SMP, magnitude, time_imaginary_part, cepstrum, frequency_imaginary_part);
        for ( m = 0; m <(ANALYSIS_WINDOW_SIZE_IN_SMP); m++) {
            fprintf(fp_save_C, "%.15f\t",cepstrum[m]);
            fflush(fp_save_C);
        }
        for (m = (CEPSTRUM_ORDER+1); m < (ANALYSIS_WINDOW_SIZE_IN_SMP-CEPSTRUM_ORDER); m++) {
            cepstrum[m] = 0 ;
        }
        n += frame_size_in_smp;
    }
    fclose(fp_save);
    fclose(fp_save_C);
    fclose(fp_save_Reespec);
    
    HammingClean();
    free(pcmbuffer);
    
    return 1;
    
    
    
    
}
