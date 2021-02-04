/*
 * Copyright 2021 Benjamin Edgington
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "c_kzg.h"

// General Utilities
void print_bytes_as_hex(byte *bytes, int start, int len);
void print_bytes_as_hex_le(byte *bytes, int start, int len);

// Fr utilities
void print_fr(const blst_fr *a);

// Fp Utilities
void print_limbs(const blst_fp *fp);

// G1 and G2 utilities
void print_p1_bytes(byte p1[96]);
void print_p1(const blst_p1 *p1);
void print_p1_affine(const blst_p1_affine *p1);
void print_p1_limbs(const blst_p1 *p1);
void print_p1_affine_limbs(const blst_p1_affine *p1);
void print_p2_affine(const blst_p2_affine *p2);
