//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>

Dna::Dna(int dna_length, Threefry::Gen &&rng)
{
    n_bytes_ = ceil((dna_length + PROM_SIZE) / 8.0); // number of words bytes to store the full DNA + mirrored beginning
    seq_ = std::vector<uint8_t>(n_bytes_);

    length_ = dna_length;

    // Generate a random genome
    size_t byte_idx;
    size_t bit_idx;

    // Fill DNA with random data
    for (int i = 0; i < dna_length; ++i)
    {
        unsigned int x = rng.random(NB_BASE);

        byte_idx = i / 8;
        bit_idx = i % 8;
        seq_[byte_idx] |= x << bit_idx;

        // Mirror bit
        if (i < PROM_SIZE)
        {
            byte_idx = (i + dna_length) / 8;
            bit_idx = (i + dna_length) % 8;
            seq_[byte_idx] |= x << bit_idx;
        }
    }
}

void Dna::print() const
{
    for (int i = 0; i < length_; ++i)
    {
        if (i % 8 == 0 && i != 0)
        {
            std::cout << " ";
        }
        std::cout << (bool)(seq_[i / 8] >> (i % 8) & uint8_t(1));
    }
    for (int i = length_; i < length_ + PROM_SIZE; i++)
    {
        if (i % 8 == 0 && i != 0)
        {
            std::cout << " ";
        }
        std::cout << (bool)(seq_[i / 8] >> (i % 8) & uint8_t(1));
    }
    std::cout << std::endl;
}

int Dna::length() const
{
    return length_;
}

void Dna::save(gzFile backup_file)
{
    size_t n_bytes = n_bytes_;
    gzwrite(backup_file, &n_bytes, sizeof(n_bytes));
    gzwrite(backup_file, seq_.data(), n_bytes * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file)
{
    size_t n_bytes;
    gzread(backup_file, &n_bytes, sizeof(n_bytes));

    uint8_t tmp_seq[n_bytes];
    gzread(backup_file, tmp_seq, n_bytes * sizeof(tmp_seq[0]));

    seq_ = std::vector<uint8_t>(tmp_seq, tmp_seq + n_bytes);
}

void Dna::set(int pos, bool val)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    // Clear bit, then set its value
    seq_[byte_idx] &= ~(1 << bit_idx);
    seq_[byte_idx] |= val << bit_idx;

    // Mirror bit in the last region of the DNA (22 lsb)
    if (pos < PROM_SIZE)
    {
        byte_idx = (pos + length_) / 8;
        bit_idx = (pos + length_) % 8;
        seq_[byte_idx] &= ~(1 << bit_idx);
        seq_[byte_idx] |= val << bit_idx;
    }
}

void Dna::do_switch(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    seq_[byte_idx] ^= 1 << bit_idx;

    // Mirror bit in the last region of the DNA (22 lsb)
    if (pos < PROM_SIZE)
    {
        byte_idx = (pos + length_) / 8;
        bit_idx = (pos + length_) % 8;
        seq_[byte_idx] ^= 1 << bit_idx;
    }
}

int Dna::promoter_at(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    // Retrieve 4 bytes of DNA including the 22 bits that interests us
    uint32_t dna_partial = *(uint32_t *)(seq_.data() + byte_idx);

    // Chop chop 🔪 the DNA to 22 bits and calculate its Hamming distance to the promoter sequence:
    // - First we aligning the bit that interests us (the bit at position `pos` in the DNA) on the right hand side of our 4 bytes.
    // - Then we only keeping the 22 least significant bits so that we can easily compute the Hamming distance to the promoter sequence.
    uint32_t dna_to_cmp = (dna_partial >> bit_idx) & PROM_MASK;

    return __builtin_popcount(dna_to_cmp ^ PROM_SEQ);
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X X !d !c !b !a
int Dna::terminator_at(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    // Retrieve 4 bytes of DNA including the 11 bits that interests us
    // We are obliged to get 4 bytes because worst case scenario is the position being on the right
    uint32_t dna_partial = *(uint32_t *)(seq_.data() + byte_idx);
    uint16_t dna_to_cmp = (dna_partial >> bit_idx) & TERM_MASK;
    return TERM_DIST_LOOKUP[dna_to_cmp];
}

bool Dna::shine_dal_start(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    uint32_t dna_partial = *(uint32_t *)(seq_.data() + byte_idx);
    uint16_t dna_to_cmp = (dna_partial >> bit_idx) & SHINE_DAL_SEQ_MASK;
    uint8_t dist = __builtin_popcount(dna_to_cmp ^ SHINE_DAL_SEQ);

    // The sequence is only considered a start if it is strictly equal to the Shine-Dalgarno sequence 011011****000
    return dist == 0;
}

bool Dna::protein_stop(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    uint16_t dna_partial = *(uint16_t *)(seq_.data() + byte_idx);
    uint16_t dna_to_cmp = (dna_partial >> bit_idx) & PROTEIN_END_MASK;
    uint8_t dist = __builtin_popcount(dna_to_cmp ^ PROTEIN_END);

    return dist == 0;
}

int Dna::codon_at(int pos)
{
    size_t byte_idx = pos / 8;
    size_t bit_idx = pos % 8;

    uint16_t dna_partial = *(uint16_t *)(seq_.data() + byte_idx);
    uint8_t codon = (dna_partial >> bit_idx) & CODON_MASK;

    return CODON_VALUE_LOOKUP[codon];
}