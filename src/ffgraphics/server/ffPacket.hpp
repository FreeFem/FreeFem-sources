#ifndef FF_PACKET_HPP
#define FF_PACKET_HPP

#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include "RNM.hpp"

using json = nlohmann::json;

template <typename T, class K>
struct ffFE {

    ffFE(KN<T> psub, KN<int> ksub, KN<K> v1)
        : Psub(psub), Ksub(ksub), V1(v1) { }

    KN_<T> Psub;
    KN_<int> Ksub;
    KN_<K> V1;
};

struct ffPacket {
    // Basic Constructor
    ffPacket();
    ffPacket(json data);

    // Copy Constructor
    ffPacket(const ffPacket& copy);
    ffPacket& operator=(const ffPacket& copy);

    // Methods
    std::string dump(int indent = -1);
    void compress();
    void clear();

    template <typename T>
    void jsonify(const T& data);

    json m_JSON;
    std::vector<uint8_t> m_Data;
    json m_Header;
};

#endif // FF_PACKET_HPP