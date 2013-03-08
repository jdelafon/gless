#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gless.h"
// gcc main.c -o main.o -c && make && ./gless


int main (int argc, const char * argv[]) {
    char track_name[] = "track_name";
    read_track(&track_name);
    printf("Hello, World!\n");
    return 0;
}


int read_track(char track_name[]){
    FILE* fichier = NULL;
    fichier = fopen(*track_name, "r+");
    if (fichier != NULL)
        {

        }
    else
        {
            printf("Impossible d'ouvrir le fichier test.txt");
        }
    fclose(fichier);

    printf("%s\n",track_name);
    return 0;
}

