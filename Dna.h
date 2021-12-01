//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <vector>
#include <zlib.h>
#include <bitset>
#include <cstdint>

#include "Threefry.h"
#include "aevol_constants.h"

class Dna
{

public:
    Dna() = default;

    Dna(const Dna &clone) = default;

    Dna(int length, Threefry::Gen &&rng);

    ~Dna() = default;

    void print() const;

    int length() const;

    void save(gzFile backup_file);

    void load(gzFile backup_file);

    void set(int pos, bool val);

    void do_switch(int pos);

    int promoter_at(int pos);

    int terminator_at(int pos);

    bool shine_dal_start(int pos);

    bool protein_stop(int pos);

    int codon_at(int pos);

private:
    size_t length_;
    size_t n_bytes_;
    std::vector<uint8_t> seq_;
};
