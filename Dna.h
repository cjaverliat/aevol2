//
// Created by arrouan on 01/10/18.
//

#pragma once

#include "Globals.h"

#include <vector>
#include <zlib.h>
#include <bitset>

#include "Threefry.h"
#include "aevol_constants.h"

class Dna {

public:
    Dna() = default;

    Dna(const Dna &clone) = default;

    Dna(int length, Threefry::Gen &&rng);

    ~Dna() = default;

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

    std::bitset<GENOME_SIZE> seq_;
};
