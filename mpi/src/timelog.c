/*
MIT License

Copyright (c) 2021 Shunji Uno <shunji_uno@iscpc.jp>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <stdio.h>
#include "timelog.h"

#ifdef _TIMELOG
void _timelog_start(timelog_t* tl) {
    gettimeofday(&(tl->stv), NULL);
}

void _timelog_end(timelog_t* tl, const char* str) {
    gettimeofday(&(tl->etv), NULL);
    float diff = (float)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((float)(tl->etv.tv_usec - tl->stv.tv_usec))/1000000.0;
    printf("TIME: %s : %f [sec]\n", str, diff);
}

float _timelog_gettime(timelog_t* tl) {
    gettimeofday(&(tl->etv), NULL);
    float diff = (float)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((float)(tl->etv.tv_usec - tl->stv.tv_usec))/1000000.0;
    return diff;
}
#endif
