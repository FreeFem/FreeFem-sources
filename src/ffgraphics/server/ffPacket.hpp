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
    ffPacket();
    ffPacket(json data);

    ffPacket(const ffPacket& copy);
    ffPacket& operator=(const ffPacket& copy);

    std::string Dump(int indent = -1);
    void Compress();
    void Clear();

    std::string GetHeader();

    template <typename T>
    void Jsonify(const T& data);

    json m_JSON;
    std::vector<uint8_t> m_Data;
    json m_Header;
};

#endif // FF_PACKET_HPP