#include "stdlib.h"
#include "stdio.h"
#include "evio.h"
#include "string.h"

#define BUFFER_LENGTH 100000

int main(int argc, char *argv[])
{
    char *ptr;
    char outputf[256], inputf[256];
    int i;
    for(i = 1; i < argc; ++i)
    {
        ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                strcpy(outputf, argv[++i]);
                break;
            case 'i':
                strcpy(inputf, argv[++i]);
                break;
            default:
                printf("Unkown option!\n");
                exit(1);
            }
        }
    }

    int inHandle, outHandle;
    int status;
    uint32_t buffer[BUFFER_LENGTH];

    status = evOpen(inputf, "r", &inHandle);
    if(status != S_SUCCESS) {
        printf("Error in open input file!\n");
        exit(1);
    }
    status = evOpen(outputf, "s", &outHandle);
    if(status != S_SUCCESS) {
        printf("Error in open output file!\n");
        exit(1);
    }

    printf("Start to split file %s to %s...\n", inputf, outputf);

    while(evRead(inHandle, buffer, BUFFER_LENGTH) == S_SUCCESS)
    {
        evWrite(outHandle, buffer);
    }

    return 0;
}
