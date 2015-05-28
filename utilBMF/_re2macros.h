/*
This extension written for re2 by Michael Axiak. Seeing if I can get speed improvements with the right c string.

Copyright (c) 2010, Michael Axiak <mike@axiak.net>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    Status API Training Shop Blog About
*/

#ifndef __RE2MACROS_H
#define __RE2MACROS_H

#include <stdio.h>
#include "re2/stringpiece.h"

static inline re2::StringPiece * new_StringPiece_array(int n)
{
    re2::StringPiece * sp = new re2::StringPiece[n];
    return sp;
}
static inline void delete_StringPiece_array(re2::StringPiece* ptr)
{
    delete[] ptr;
}

#define addressof(A) (&A)
#define addressofs(A) (&A)

#define as_char(A) (char *)(A)
#define pattern_Replace(A, B, C) re2::RE2::Replace((A), (B), (C))
#define pattern_GlobalReplace(A, B, C) re2::RE2::GlobalReplace((A), (B), (C))

#endif
