#include <iostream>
#include <vector>

enum Tag {CONJUNCTION = 0, NOUN = 1, VERB = 2, ADJECTIVE = 3};
enum Word {THAT = 0, IS = 1, NOT = 2, IT = 3, GOOD = 4};

const float probVector[] = 
    {1.f/7.f, 1.f/3.f, 3.f/7.f, 2.f/21.f};

const float transMat[][4] = {
    {1.f/6.f, 1.f/2.f, 1.f/6.f, 1.f/6.f},
    {1.f/9.f, 2.f/9.f, 4.f/9.f, 2.f/9.f},
    {1.f/10.f, 4.f/10.f, 4.f/10.f, 1.f/10.f},
    {1.f/4.f, 1.f/4.f, 1.f/4.f, 1.f/4.f}};

const float emissionMat[][5] = {
    {3.f/7.f, 1.f/7.f, 1.f/7.f, 1.f/7.f, 1.f/7.f},
    {4.f/11.f, 1.f/11.f, 2.f/11.f, 3.f/11.f, 1.f/11.f},
    {2.f/13.f, 7.f/13.f, 2.f/13.f, 1.f/13.f, 1.f/13.f},
    {1.f/6.f, 1.f/6.f, 1.f/6.f, 1.f/6.f, 2.f/6.f}};

// Forward algorithm
float observed_seq_prob(const std::vector<Word>& observed_seq,
        std::vector<std::vector<float>>& prob_vec, bool verbose = false) {
    prob_vec.push_back({0.f, 0.f, 0.f, 0.f});
    int dyn_index = 0; // dynamic index
    
    // initial computation
    for (int i = 0; i < 4; i++) {
        prob_vec[dyn_index][i] = probVector[i] * emissionMat[i][observed_seq[dyn_index]];
        if (verbose) std::cout << "prob_vec[" << dyn_index << "][" << i << "] = " << prob_vec[dyn_index][i] << '\n';
    }
    // rest of the sequence, using previously memoized results
    for (int i = 1; i < observed_seq.size(); i++) {
        dyn_index++;
        prob_vec.push_back({0.f, 0.f, 0.f, 0.f});
        for (int j = 0; j < 4; j++) {
            prob_vec[dyn_index][j] = prob_vec[dyn_index-1][j] * emissionMat[j][i];
            for (int k = 0; k < 4; k++) {
                if (j == k) continue; // ignore self-transition
                prob_vec[dyn_index][j] += transMat[k][j] * prob_vec[dyn_index-1][k];
            }
            prob_vec[dyn_index][j] *= emissionMat[j][observed_seq[i]];
            if (verbose) std::cout << "prob_vec[" << dyn_index << "][" << j << "] = " << prob_vec[dyn_index][j] << '\n';
        }
    }
    
    float ret = 0;
    for (const auto& p : prob_vec[dyn_index]) {
        ret += p;
    }
    return ret;
}

float prob_of_tag_given_seq(Tag tag, float seq_prob, std::vector<std::vector<float>>& prob_vec, bool verbose = false) {
    float prob = 0.f;
    int index = prob_vec.size() - 1;
    
    for (int i = 0; i < 4; i++) {
        prob += transMat[i][tag] * (prob_vec[index][tag] / seq_prob);
    }
    
    return prob;
}

// Viterbi's Algorithm
std::vector<Tag> most_likely_tags(const std::vector<Word>& observed_seq, std::vector<std::vector<float>>& prob_vec, bool verbose = false) {
    prob_vec.push_back({0.f, 0.f, 0.f, 0.f});
    int dyn_index = 0; // dynamic index
    
    // initial computation
    for (int i = 0; i < 4; i++) {
        prob_vec[dyn_index][i] = probVector[i] * emissionMat[i][observed_seq[dyn_index]];
        if (verbose) std::cout << "prob_vec[" << dyn_index << "][" << i << "] = " << prob_vec[dyn_index][i] << '\n';
    }
    // rest of the sequence, using previously memoized results
    for (int i = 1; i < observed_seq.size(); i++) {
        dyn_index++;
        prob_vec.push_back({0.f, 0.f, 0.f, 0.f});
        for (int j = 0; j < 4; j++) {
            auto max_p = prob_vec[dyn_index-1][j] * emissionMat[j][i];
            for (int k = 0; k < 4; k++) {
                if (j == k) continue; // ignore self-transition
                auto t = transMat[k][j] * prob_vec[dyn_index-1][k];
                if (t > max_p) max_p = t;
            }
            prob_vec[dyn_index][j] = max_p * emissionMat[j][observed_seq[i]];
            if (verbose) std::cout << "prob_vec[" << dyn_index << "][" << j << "] = " << prob_vec[dyn_index][j] << '\n';
        }
    }
    
    std::vector<Tag> tags;
    for (int i = 0; i < prob_vec.size(); i++) {
        int max_i = 0;
        for (int j = 1; j < prob_vec[0].size(); j++) {
            if (prob_vec[i][j] > prob_vec[i][max_i]) {
                max_i = j;
            }
        }
        tags.push_back(static_cast<Tag>(max_i));
    }
    return tags;
}

std::string tags_to_string(std::vector<Tag> tags) {
    std::string ret = "";
    for (const auto& t : tags) {
        switch (t) {
            case CONJUNCTION:   ret += " Conjunction"; break;
            case NOUN:          ret += " Noun"; break;
            case VERB:          ret += " Verb"; break;
            case ADJECTIVE:     ret += " Adjective"; break;
        }
    }
    return ret;
}

std::string words_to_string(std::vector<Word> words) {
    std::string ret = "";
    for (const auto& w : words) {
        switch (w) {
            case THAT:  ret += " That"; break;
            case IS:    ret += " Is"; break;
            case NOT:   ret += " Not"; break;
            case IT:    ret += " It"; break;
            case GOOD:  ret += " Good"; break;
        }
    }
    return ret;
}

int main(int argc, char* argv[]) {
    bool v = false;
    if (argc > 1) {
        if (std::string(argv[1]) == "-v") v = true;
    }
    
    std::vector<std::vector<float>> prob_vec; // dynamic 2D array to hold computation results between steps
    std::vector<Word> seq_not_that_good;
    seq_not_that_good.push_back(NOT);
    seq_not_that_good.push_back(THAT);
    seq_not_that_good.push_back(GOOD);
    std::cout << "---------- QUESTION 1 b) ----------" << '\n';
    auto prob = observed_seq_prob(seq_not_that_good, prob_vec, v);
    std::cout << "The probability for" << words_to_string(seq_not_that_good) << " is: " << prob << '\n';
    std::cout << '\n' << "---------- QUESTION 1 c) ----------" << '\n';
    auto prob2 = prob_of_tag_given_seq(NOUN, prob, prob_vec, v);
    std::cout << "The probability of having a noun after" << words_to_string(seq_not_that_good) << " is: " << prob2 << '\n';
    std::cout << '\n' << "---------- QUESTION 1 d) ----------" << '\n';
    prob_vec.clear();
    auto tags = most_likely_tags(seq_not_that_good, prob_vec, v);
    std::cout << "The most likey tags for" << words_to_string(seq_not_that_good) << " is:" << tags_to_string(tags) << '\n';
    
    return 0;
}

