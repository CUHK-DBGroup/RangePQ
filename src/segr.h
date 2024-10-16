#ifndef SEGR_H
#define SEGR_H

#include <iostream>
#include <cassert>
#include "pqkmeans.h"
#include "distance.h"
#include "WBTTree.hpp"
// #include "WBTChunk.hpp"
#include "PQRS.hpp"
#include "BST.hpp"
#include <unordered_set>
#include <unordered_map>
#include <sys/resource.h>
#include <fstream>
#include <map>
#include <chrono>

// Handle missing ssize_t on Windows. 
# if defined(_MSC_VER) 
    typedef __int64 ssize_t;
# endif

// For py::array_t
// See http://pybind11.readthedocs.io/en/master/advanced/pycpp/numpy.html#direct-access
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;


namespace segr {

struct DistanceTable{
    // Helper structure. This is identical to vec<vec<float>> dt(M, vec<float>(Ks))
    DistanceTable() {}
    DistanceTable(size_t M, size_t Ks) : Ks_(Ks), data_(M * Ks) {}
    void SetVal(size_t m, size_t ks, float val) {
        data_[m * Ks_ + ks] = val;
    }
    float GetVal(size_t m, size_t ks) const {
        return data_[m * Ks_ + ks];
    }
    size_t Ks_;
    std::vector<float> data_;
};




class SegrCpp {
public:
    std::vector<int> coarse_id;

    // WBTRotate *WBT = nullptr;
    WBTChunk *WBT = nullptr;

    int chunkSize = 2000;
    SegrCpp() {}  // Shouldn't be default-constructed
    SegrCpp(const py::array_t<float> &codewords, bool verbose);

    // ===== Functions that can be called from Python =====
    //void SetCodewords(const py::array_t<float> &codewords);  // This should be called first
    void Reconfigure(int nlist, int iter, int cz);
    void Reconfigure_2(int nlist, int iter, int cz, const py::array_t<unsigned char> &codes,
                    const py::array_t<size_t> &ass);
    void AddCodes(const py::array_t<unsigned char> &codes,
                    const py::array_t<size_t> &keySet,
                    const py::array_t<int> &valueSet,
                    bool update_flag);
    // The default integers of Python is int64 (long long), so the type of target_ids is long long
    void DeleteKey(size_t key);
    std::pair<std::vector<size_t>, std::vector<float>> QueryLinear(const py::array_t<float> &query,
                                                                int topk,
                                                                const py::array_t<long long> &target_ids) const;
    void Clear();

    // ===== Functions that would not be called from Python (Used inside c++) =====
    void UpdatePostingLists(size_t start, size_t num);
    void UpdatePostingLists_2(size_t start, size_t num, const py::array_t<size_t> &ass);
    DistanceTable DTable(const py::array_t<float> &vec) const;
    float ADist(const DistanceTable &dtable, const std::vector<unsigned char> &code) const;
    float ADist(const DistanceTable &dtable, const std::vector<unsigned char> &flattened_codes, size_t n) const;
    std::pair<std::vector<size_t>, std::vector<float>> PairVectorToVectorPair(const std::vector<std::pair<size_t, float>> &pair_vec) const;
    std::pair<std::vector<size_t>, std::vector<float> > QuerySegIVF(const py::array_t<float>& query,
        int topk,
        std::unordered_set<int>& range_union_set,
        int range_l, int range_r,
        int L)const;
    // Property getter
    size_t GetN() const {return flattened_codes_.size() / M_;}
    size_t GetNumList() const {return coarse_centers_.size();}

    // Given a long (N * M) codes, pick up n-th code
    std::vector<unsigned char> NthCode(const std::vector<unsigned char> &long_code, size_t n) const;
    // Given a long (N * M) codes, pick up m-th element from n-th code
    unsigned char NthCodeMthElement(const std::vector<unsigned char> &long_code, std::size_t n, size_t m) const;

    std::pair<std::vector<size_t>, std::vector<float> > QuerySegRange(const py::array_t<float>& query,
        int topk,
        int range_l, int range_r,
        int L);
    // Member variables
    // M_ -> the number of new dimensions
    // Ks_ -> the number of identifier
    size_t M_, Ks_;
    bool verbose_;

    // M-dimensions Ks pre-trained code words (identifier) Ds sub-vector
    std::vector<std::vector<std::vector<float>>> codewords_;  // (M, Ks, Ds)
    std::vector<std::vector<unsigned char>> coarse_centers_;  // (NumList, M)
    std::vector<unsigned char> flattened_codes_;  // (N, M) PQ codes are flattened to N * M long array

    std::vector<BST* > postingList;

    std::vector<int> keyList;
    std::vector<int> valueList;
    std::unordered_map<int, int> key2Id;
    std::vector<std::pair<std::pair<int,float>, int > > w;

    float query_range_selectivity(int range_l, int range_r);

    void Memory(){

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        float gms = ms / 1024 ;
        std::cout <<"Memory " << gms << "mb\n";

        std::ifstream file("/proc/self/status");
        std::string line;

        while (std::getline(file, line)) {
            if (line.find("VmRSS") != std::string::npos) {
                std::cout << line << std::endl;
                break;
            }
        }
    }
 
    WBTChunk *buildWBT(){
        return new WBTChunk(valueList, keyList, coarse_id, chunkSize);
    }

    void updateWBT(size_t start, size_t num){
        for(size_t n = 0; n < num; n++){
            WBT->insert(valueList[start + n], keyList[start + n], coarse_id[start + n]);
        }
    }
};


SegrCpp::SegrCpp(const py::array_t<float> &codewords, bool verbose)
{
    verbose_ = verbose;
    const auto &r = codewords.unchecked<3>();  // codewords must have ndim=3, with non-writable
    M_ = (size_t) r.shape(0);
    Ks_ = (size_t) r.shape(1);
    size_t Ds = (size_t) r.shape(2);
    codewords_.resize(M_, std::vector<std::vector<float>>(Ks_, std::vector<float>(Ds)));
    for (ssize_t m = 0; m < r.shape(0); ++m) {
        for (ssize_t ks = 0; ks < r.shape(1); ++ks) {
            for (ssize_t ds = 0; ds < r.shape(2); ++ds) {
                codewords_[m][ks][ds] = r(m, ks, ds);
            }
        }
    }

    // if (verbose_) {
    //     // Check which SIMD functions are used. See distance.h for this global variable.
    //     std::cout << "SIMD support: " << g_simd_architecture << std::endl;
    // }
}
// nlist posting list number [k in paper]
void SegrCpp::Reconfigure(int nlist, int iter, int cz)
{
    chunkSize = cz;

    assert(0 < nlist);
    assert((size_t) nlist <= GetN());

    // ===== (1) Sampling vectors for pqk-means =====
    // Since clustering takes time, we use a subset of all codes for clustering.
    size_t len_for_clustering = std::min(GetN(), (size_t) nlist * 100);
    // if (verbose_) {
    //     std::cout << "The number of vectors used for training of coarse centers: " << len_for_clustering << std::endl;
    // }
    // Prepare a random set of integers, drawn from [0, ..., N-1], where the cardinality of the set is len_for_clustering
    std::vector<size_t> ids_for_clustering(GetN());  // This can be large and might be the bootle neck of memory consumption
    std::iota(ids_for_clustering.begin(), ids_for_clustering.end(), 0);  // 0, 1, 2, ...
    std::shuffle(ids_for_clustering.begin(), ids_for_clustering.end(), std::default_random_engine(123));
    ids_for_clustering.resize(len_for_clustering);
    ids_for_clustering.shrink_to_fit();  // For efficient memory usage

    std::vector<unsigned char> flattened_codes_randomly_picked;  // size=len_for_clustering
    flattened_codes_randomly_picked.reserve(len_for_clustering * M_);
    for (const auto &id : ids_for_clustering) {  // Pick up vectors to construct a training set
        std::vector<unsigned char> code = NthCode(flattened_codes_, id);
        flattened_codes_randomly_picked.insert(flattened_codes_randomly_picked.end(),
                                               code.begin(), code.end());
    }
    assert(flattened_codes_randomly_picked.size() == len_for_clustering * M_);


    // ===== (2) Run pqk-means =====
    if (verbose_) {std::cout << "Start to run PQk-means" << std::endl;}
    pqkmeans::PQKMeans clustering_instance(codewords_, nlist, iter, verbose_);
    clustering_instance.fit(flattened_codes_randomly_picked);


    // ===== (3) Update coarse centers =====
    coarse_centers_ = clustering_instance.GetClusterCenters();
    assert(coarse_centers_.size() == (size_t) nlist);
    assert(coarse_centers_[0].size() == M_);


    // ===== (4) Update posting lists =====
    if (verbose_) {std::cout << "Start to update posting lists" << std::endl;}

    // for(int i = 0; i < postingList.size(); i++)
    //     if(postingList[i] != nullptr)
    //         postingList[i]->clear();
    // postingList.clear();
    // postingList.resize(nlist);
    // for(int i = 0; i < nlist; i++)
    //     postingList[i] = new BST();

    
        struct rusage rUsage2;
        getrusage(RUSAGE_SELF, &rUsage2);
        long ms2 = rUsage2.ru_maxrss;
        float gms2 = ms2 / 1024.0 ;
        std::cout <<"Memory " << gms2 << "mb\n";
    UpdatePostingLists(0, GetN());

    puts("updateOK");

    puts("buildBegin");
    if(WBT != nullptr) delete WBT;

        // struct rusage rUsage;
        // getrusage(RUSAGE_SELF, &rUsage);
        // long ms = rUsage.ru_maxrss;
        // float gms = ms / 1024.0 ;
        // std::cout <<"Memory " << gms << "mb\n";
        std::ifstream file("/proc/self/status");
        std::string line;

        while (std::getline(file, line)) {
            if (line.find("VmRSS") != std::string::npos) {
                std::cout << line << std::endl;
                break;
            }
        }

    WBT =buildWBT();
    std::cout<<WBT->chunkSize<<std::endl;
    std::unordered_set<int> range_union_set;
    
    std::ifstream file1("/proc/self/status");
        std::string line1;
        while (std::getline(file1, line1)) {
            if (line1.find("VmRSS") != std::string::npos) {
                std::cout << line1 << std::endl;
                break;
            }
        }
        // struct rusage rUsage1;
        // getrusage(RUSAGE_SELF, &rUsage1);
        // long ms1 = rUsage1.ru_maxrss;
        // float gms1 = ms1/ 1024.0 ;
        // std::cout <<"Memory " << gms1 - gms << "mb\n";


    puts("buildOK");

    // struct rusage sage;
    // getrusage(RUSAGE_SELF, &sage);
    // ms = sage.ru_maxrss;
    // gms = ms / 1024 ;
    // std::cout <<"memory " << gms << "mb\n";

}


void SegrCpp::DeleteKey(size_t key){
    int id = key2Id[key];
    // auto t1 = std::chrono::steady_clock::now();
    // printf("%d %d", keyList[id], valueList[id]);
    WBT->erase(valueList[id], keyList[id], coarse_id[id]);
    // printf("erase %ld",key);puts("OK");
    // auto t2 = std::chrono::steady_clock::now();
    // // postingList[coarse_id[id]]->erase(keyList[id], valueList[id]);
    // auto t3 = std::chrono::steady_clock::now();
    int N = keyList.size();
    key2Id[keyList[N-1]] = id;
    key2Id.erase(key);
    std::swap(keyList[id], keyList[N - 1]);
    keyList.pop_back();
    std::swap(valueList[id], valueList[N - 1]);
    valueList.pop_back();
    std::swap(coarse_id[id], coarse_id[N - 1]);
    coarse_id.pop_back();
    for(int i = 0; i < M_; i++){
        std::swap(flattened_codes_[ (id) * M_ + i], flattened_codes_[ (N - 1) * M_ + i]);
    }
    for(int i = 0; i < M_; i++){
        flattened_codes_.pop_back();
    }
    // printf("erase %ld",key);puts("OK2");
    // auto t4 = std::chrono::steady_clock::now();
    // auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1);
    // auto d2 = std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2);
    // auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3);
    // std::cout<<d1.count()<<' '<<d2.count()<<' '<<d3.count()<<std::endl;
}

void SegrCpp::AddCodes(const py::array_t<unsigned char> &codes,
                    const py::array_t<size_t> &keySet,
                    const py::array_t<int> &valueSet,
                    bool update_flag)
{
    // (1) Add new input codes to flatted_codes. This imply pushes back the elements.
    // After that, if update_flg=true, (2) update posting lists for the input codes.
    // Note that update_flag should be true in usual cases. It should be false
    // if (1) this is the first call of AddCodes (i.e., calling in add_configure()),
    // of (2) you've decided to call reconfigure() manually after add()

    if (update_flag && coarse_centers_.empty()) {
        std::cerr << "Error. reconfigure() must be called before running add(vecs=X, update_posting_lists=True)."
                  << "If this is the first addition, please call add_configure(vecs=X)" << std::endl;
        throw;
    }

    // ===== (1) Add codes to flattened_codes =====
    const auto &r = codes.unchecked<2>(); // codes must have ndim=2; with non-writeable
    const auto &k = keySet.unchecked<1>();
    const auto &v = valueSet.unchecked<1>();
    size_t N = (size_t) r.shape(0);
    assert(M_ == (size_t) r.shape(1));
    size_t N0 = GetN();
    flattened_codes_.resize( (N0 + N) * M_);

    keyList.resize(N0 + N);
    valueList.resize(N0 + N);

    for (size_t n = 0; n < N; ++n) {
        keyList[N0 + n] = k(n);
        valueList[N0 + n] = v(n);
        key2Id[k(n)] = N0 + n;
        for (size_t m = 0; m < M_; ++m) {
            flattened_codes_[ (N0 + n) * M_ + m] = r(n, m);
        }
    }
    if (verbose_) {
        std::cout << N << " new vectors are added." << std::endl;
        std::cout << "Total number of codes is " << GetN() << std::endl;
    }

    // ===== (2) Update posting lists =====
    if (update_flag) {
        if (verbose_) { std::cout << "Start to update posting lists" << std::endl; }
        UpdatePostingLists(N0, N);
        updateWBT(N0, N);
    }
}

std::pair<std::vector<size_t>, std::vector<float> > SegrCpp::QueryLinear(const py::array_t<float> &query,
                                                                        int topk,
                                                                        const py::array_t<long long> &target_ids) const
{
    const auto &tids = target_ids.unchecked<1>(); // target_ids must have ndim = 1; can be non-writeable
    size_t S = tids.shape(0);  // The number of target_ids. It might be 0 if not specified.

    assert((size_t) topk <= GetN());

    // ===== (1) Create dtable =====
    DistanceTable dtable = DTable(query);

    // ===== (2) Run PQ linear search =====
    // [todo] Can be SIMDized?
    std::vector<std::pair<size_t, float>> scores;
    if (S == 0) {  // No target ids
        size_t N = GetN();
        scores.resize(N);
// #pragma omp parallel for
        for (long long n_tmp = 0LL; n_tmp < static_cast<long long>(N); ++n_tmp) {
            size_t n = static_cast<size_t>(n_tmp);
            scores[n] = {n, ADist(dtable, flattened_codes_, n)};
        }
    } else {  // Target ids are specified
        assert((size_t) topk <= S);
        assert(S <= GetN());
        scores.resize(S);
// #pragma omp parallel for
        for (long long s_tmp = 0LL; s_tmp < static_cast<long long>(S); ++s_tmp) {
            size_t s = static_cast<size_t>(s_tmp);
            size_t tid = static_cast<size_t>(tids(s));
            scores[s] = {tid, ADist(dtable, flattened_codes_, tid)};
        }
    }



    // ===== (3) Sort them =====
    // [todo] Can be parallelized?
    std::partial_sort(scores.begin(), scores.begin() + topk, scores.end(),
                      [](const std::pair<size_t, float> &a, const std::pair<size_t, float> &b){return a.second < b.second;});
    scores.resize(topk);
    scores.shrink_to_fit();

    // ===== (4) Return the result, in the form of pair<vec, vec> =====
    // Note that this returns two lists, not np.array
    return PairVectorToVectorPair(scores);
}


float SegrCpp::query_range_selectivity(int range_l, int range_r){
    int test_num = GetN() * 0.001;
    int tot_success_num = 0;
    int totNum = GetN();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, totNum-1);


    for(int i = 1;i<=test_num;i++){
        
        int random_number = distrib(gen);
        if(valueList[random_number]<=range_r&&valueList[random_number]>=range_l){
            tot_success_num ++;
        }
    }
    // std::cout<<"Query_Selectivity"<< " "<<range_l<<" "<<range_r<<" " <<(float) tot_success_num / test_num<<"\n";
    return (float) tot_success_num / test_num;
}

std::pair<std::vector<size_t>, std::vector<float> > SegrCpp::QuerySegRange(const py::array_t<float>& query,
    int topk,
    int range_l, int range_r,
    int L) {
    std::unordered_set<int> range_union_set;
    
    WBT->ask(range_union_set, range_l, range_r);
    // printf("%d %d",range_l,range_r);
    // puts("WBTOK");
    // for(auto x:range_union_set)
    //     printf("%d ",x);
    // std::cout << "Hello, world!" << std::flush;

    // printf("%d\n",L);
    L = std::max((int)(L * query_range_selectivity(range_l, range_r)), 1000);
    // printf("%d\n",L);
    auto res = QuerySegIVF(query, topk, range_union_set, range_l, range_r, L);
    // auto t3 = std::chrono::steady_clock::now();

    // auto d1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    // auto d2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2);

    // std::cout << d1.count() << " " << d2.count()<<" "<<std::endl;
    return res;
}

std::pair<std::vector<size_t>, std::vector<float> > SegrCpp::QuerySegIVF(const py::array_t<float> &query,
                                                                     int topk,
                                                                     std::unordered_set<int> &range_union_set,
                                                                     int range_l, int range_r,
                                                                     int L) const {
    //const auto &tids = target_ids.unchecked<1>(); // target_ids must have ndim = 1 with non-writeable
    //size_t S = tids.shape(0);  // The number of target_ids. It might be 0 if not specified.

    assert((size_t) topk <= GetN());
    assert(topk <= L && (size_t) L <= GetN());
    

    // ===== (1) Create dtable =====
    DistanceTable dtable = DTable(query);      

    // ===== (2) Compare to coarse centers and sort the results =====
    std::vector<std::pair<size_t, float>> scores_coarse;
    size_t nlist = GetNumList();
    //#pragma omp parallel for
    for (auto no : range_union_set) {
        // printf("%d", no);
        // puts("");
        scores_coarse.push_back({no, ADist(dtable, coarse_centers_[no]) });
    }

    // for (size_t no = 0; no < nlist; ++no) {
    //     scores_coarse.push_back( {no, ADist(dtable, coarse_centers_[no])});
    // }

    std::sort(scores_coarse.begin(), scores_coarse.end(),
        [](const std::pair<size_t, float>& a, const std::pair<size_t, float>& b) {return a.second < b.second; });

    // puts("sort OK");

    std::vector<int> labelId;
    for(int i = 0; i < scores_coarse.size();i++) {
        // printf("%d ", scores_coarse[i].first);
        labelId.push_back(scores_coarse[i].first);
    }
    // std::cout << "Hello, world!" << std::flush;
    std::vector<int> vecId;

    // puts("");

    WBT->getId(labelId,vecId,range_l,range_r,L);

    L = std::min(L, (int)vecId.size());

    std::vector<std::pair<size_t, float>> scores;
    for(int i = 0;i<vecId.size();i++){
        if(key2Id.count(vecId[i]) == 0) continue;
        int id = key2Id.at(vecId[i]);
        // printf("%d:%d ",vecId[i], id);
        scores.emplace_back(vecId[i], ADist(dtable, flattened_codes_, id));

    }
    
    // std::cout << "Hello, world!" << std::flush;
    // puts("OK");

    if (true) {
        // ===== (8) Sort them =====
        std::partial_sort(scores.begin(), scores.begin() + topk, scores.end(),
            [](const std::pair<size_t, float>& a, const std::pair<size_t, float>& b) {return a.second < b.second; });
        scores.resize(topk);
        scores.shrink_to_fit();

        // ===== (9) Return the result, in the form of pair<vec, vec> =====
        // Note that this returns two lists, not np.array
        return PairVectorToVectorPair(scores);
    }

}

void SegrCpp::Clear()
{
    coarse_centers_.clear();
    flattened_codes_.clear();
    // for(auto p: postingList) delete p;
    // postingList.clear();
}

void SegrCpp::UpdatePostingLists(size_t start, size_t num)
{
    // Update (add) identifiers to posting lists, from codes[start] to codes[start + num -1]
    // This just add IDs, so be careful to call this (e.g., the same IDs will be added if you call
    // this funcs twice at the same time, that would be not expected behavior)
    assert(start <= GetN());
    assert(start + num <= GetN());


    // ===== (1) Construct a dummy pqkmeans class for computing Symmetric Distance =====
    pqkmeans::PQKMeans clustering_instance(codewords_, (int)GetNumList(), 0, true);
    clustering_instance.SetClusterCenters(coarse_centers_);
    std::vector<size_t> assign(num);
    coarse_id.resize(start + num);
    if(num >100){
#pragma omp parallel for
        for (long long n_tmp = 0LL; n_tmp < static_cast<long long>(num); ++n_tmp) {
            size_t n = static_cast<size_t>(n_tmp);
            assign[n] = clustering_instance.predict_one(NthCode(flattened_codes_, start + n));
        }
    }
    else{
        for (long long n_tmp = 0LL; n_tmp < static_cast<long long>(num); ++n_tmp) {
            size_t n = static_cast<size_t>(n_tmp);
            assign[n] = clustering_instance.predict_one_single(NthCode(flattened_codes_, start + n));
        }
    }

    for (size_t n = 0; n < num; ++n) {
        int id = start + n;
        // postingList[assign[n]]->insert(keyList[id], valueList[id]);
        coarse_id[start + n] = assign[n];
    }
    // if(verbose_){
    //     std::cout<<"cluster of"<<assign[1062]<<"\n";
    // }
}

void SegrCpp::UpdatePostingLists_2(size_t start, size_t num, const py::array_t<size_t> &ass)
{
    // Update (add) identifiers to posting lists, from codes[start] to codes[start + num -1]
    // This just add IDs, so be careful to call this (e.g., the same IDs will be added if you call
    // this funcs twice at the same time, that would be not expected behavior)
    assert(start <= GetN());
    assert(start + num <= GetN());

    coarse_id.resize(start + num);

    // ===== (1) Construct a dummy pqkmeans class for computing Symmetric Distance =====
    const auto &k = ass.unchecked<1>();
    for (size_t n = 0; n < num; ++n) {
        int id = start + n;
        // postingList[assign[n]]->insert(keyList[id], valueList[id]);
        coarse_id[start + n] = k(n);
    }
    // if(verbose_){
    //     std::cout<<"cluster of"<<assign[1062]<<"\n";
    // }
}

void SegrCpp::Reconfigure_2(int nlist, int iter, int cz, const py::array_t<unsigned char> &codes,
                    const py::array_t<size_t> &ass)
{
    chunkSize = cz;
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024 ;
    std::cout <<"Memory " << gms << "mb\n";

    assert(0 < nlist);
    assert((size_t) nlist <= GetN());

    // ===== (3) Update coarse centers =====
    coarse_centers_.resize(nlist);
    const auto &r = codes.unchecked<2>(); 
    for(int i = 0; i< nlist;i++)
        for(int j = 0; j < M_; j++){
            coarse_centers_[i].push_back(r(i,j));
        }
    assert(coarse_centers_.size() == (size_t) nlist);
    assert(coarse_centers_[0].size() == M_);


    // ===== (4) Update posting lists =====
    if (verbose_) {std::cout << "Start to update posting lists" << std::endl;}

    UpdatePostingLists_2(0, GetN(), ass);

    puts("updateOK");

    puts("buildBegin");
    if(WBT != nullptr) delete WBT;
    WBT =buildWBT();


    puts("buildOK");

    struct rusage sage;
    getrusage(RUSAGE_SELF, &sage);
    ms = sage.ru_maxrss;
    gms = ms / 1024 ;
    std::cout <<"memory " << gms << "mb\n";

}

DistanceTable SegrCpp::DTable(const py::array_t<float> &vec) const
{
    const auto &v = vec.unchecked<1>();
    size_t Ds = codewords_[0][0].size();
    assert((size_t) v.shape(0) == M_ * Ds);
    DistanceTable dtable(M_, Ks_);
    for (size_t m = 0; m < M_; ++m) {
        for (size_t ks = 0; ks < Ks_; ++ks) {
            dtable.SetVal(m, ks, fvec_L2sqr(&(v(m * Ds)), codewords_[m][ks].data(), Ds));
        }
    }
    return dtable;
}

float SegrCpp::ADist(const DistanceTable &dtable, const std::vector<unsigned char> &code) const
{  
    assert(code.size() == M_);
    float dist = 0;
    for (size_t m = 0; m < M_; ++m) {
        unsigned char ks = code[m];
        dist += dtable.GetVal(m, ks);
    }
    return dist;
}

float SegrCpp::ADist(const DistanceTable &dtable, const std::vector<unsigned char> &flattened_codes, size_t n) const
{
    float dist = 0;
    for (size_t m = 0; m < M_; ++m) {
        unsigned char ks = NthCodeMthElement(flattened_codes, n, m);
        dist += dtable.GetVal(m, ks);
    }
    return dist;
}

std::pair<std::vector<size_t>, std::vector<float> > SegrCpp::PairVectorToVectorPair(const std::vector<std::pair<size_t, float> > &pair_vec) const
{
    std::pair<std::vector<size_t>, std::vector<float>> vec_pair(std::vector<size_t>(pair_vec.size()), std::vector<float>(pair_vec.size()));
    for(size_t n = 0, N = pair_vec.size(); n < N; ++n) {
        vec_pair.first[n] = pair_vec[n].first;
        vec_pair.second[n] = pair_vec[n].second;
    }
    return vec_pair;
}



std::vector<unsigned char> SegrCpp::NthCode(const std::vector<unsigned char> &long_code, size_t n) const
{
    return std::vector<unsigned char>(long_code.begin() + n * M_, long_code.begin() + (n + 1) * M_);
}

unsigned char SegrCpp::NthCodeMthElement(const std::vector<unsigned char> &long_code, std::size_t n, size_t m) const
{
    return long_code[ n * M_ + m];
}


} // namespace rii

#endif // RII_H
