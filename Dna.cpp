//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"

#include <cassert>
#include <iostream>

Dna::Dna(int length, Threefry::Gen &&rng) {
    // Generate a random genome
    for (int32_t i = 0; i < length; i++) {
        seq_[length - 1 - i] = rng.random(NB_BASE);
    }
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    int dna_length = length();
    std::string seq_str = seq_.to_string();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, seq_str.data(), dna_length * sizeof(char));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(char));
    seq_ = std::bitset<GENOME_SIZE>(tmp_seq);
}

void Dna::set(int pos, bool val) {
    seq_[length() - 1 - pos] = val;
}

void Dna::do_switch(int pos) {
    seq_.flip(length() - 1 - pos);
}

int Dna::promoter_at(int pos) {
    int prom_dist[PROM_SIZE];

    for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
        int search_pos = pos + motif_id;
        if (search_pos >= seq_.size())
            search_pos -= seq_.size();
        // Searching for the promoter
        prom_dist[motif_id] =
                PROM_SEQ[PROM_SIZE - 1 - motif_id] ^ seq_[length() - 1 - search_pos];

    }


    // Computing if a promoter exists at that position
    int dist_lead = prom_dist[0] +
                    prom_dist[1] +
                    prom_dist[2] +
                    prom_dist[3] +
                    prom_dist[4] +
                    prom_dist[5] +
                    prom_dist[6] +
                    prom_dist[7] +
                    prom_dist[8] +
                    prom_dist[9] +
                    prom_dist[10] +
                    prom_dist[11] +
                    prom_dist[12] +
                    prom_dist[13] +
                    prom_dist[14] +
                    prom_dist[15] +
                    prom_dist[16] +
                    prom_dist[17] +
                    prom_dist[18] +
                    prom_dist[19] +
                    prom_dist[20] +
                    prom_dist[21];

    return dist_lead;
}

// Given a, b, c, d boolean variable and X random boolean variable,
// a terminator look like : a b c d X X !d !c !b !a

// TODO: faire pareil que promoter_at
int Dna::terminator_at(int pos) {
    int term_dist[TERM_STEM_SIZE];
    for (int motif_id = 0; motif_id < TERM_STEM_SIZE; motif_id++) {
        int right = pos + motif_id;
        int left = pos + (TERM_SIZE - 1) - motif_id;

        // loop back the dna inf needed
        if (right >= length()) right -= length();
        if (left >= length()) left -= length();

        // Search for the terminators
        term_dist[motif_id] = seq_[length() - 1 - right] ^ seq_[length() - 1 - left];
    }
    int dist_term_lead = term_dist[0] +
                         term_dist[1] +
                         term_dist[2] +
                         term_dist[3];

    return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
    bool start = false;
    int t_pos, k_t;

    for (int k = 0; k < SHINE_DAL_SIZE + CODON_SIZE; k++) {
        k_t = k >= SHINE_DAL_SIZE ? k + SD_START_SPACER : k;
        t_pos = pos + k_t;
        if (t_pos >= seq_.size())
            t_pos -= seq_.size();

        if(SHINE_DAL_SEQ[k_t] == '*') {
            start = false;
            break;
        } else if (seq_[length() - 1 - t_pos] == (SHINE_DAL_SEQ[k_t] - '0')) {
            start = true;
        } else {
            start = false;
            break;
        }
    }

    return start;
}

bool Dna::protein_stop(int pos) {
    bool is_protein;
    int t_k;

    for (int k = 0; k < CODON_SIZE; k++) {
        t_k = pos + k;
        if (t_k >= seq_.size())
            t_k -= seq_.size();

        if (seq_[length() - 1 - t_k] == (PROTEIN_END[k] - '0')) {
            is_protein = true;
        } else {
            is_protein = false;
            break;
        }
    }

    return is_protein;
}

int Dna::codon_at(int pos) {
    int value = 0;

    int t_pos;

    for (int i = 0; i < CODON_SIZE; i++) {
        t_pos = pos + i;
        if (t_pos >= seq_.size())
            t_pos -= seq_.size();
        if (seq_[length() - 1 - t_pos] == 1)
            value += 1 << (CODON_SIZE - i - 1);
    }

    return value;
}