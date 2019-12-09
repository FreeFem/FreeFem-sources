#ifndef FF_PACKET_HPP
#define FF_PACKET_HPP

#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include "RNM.hpp"

using json = nlohmann::json;

#define MAX_PACKET_SIZE 16384

template <typename T, class K>
struct ffFE {

    ffFE(KN<T> psub, KN<int> ksub, KN<K> v1, K max, K min, bool AsVector)
        : Psub(psub), Ksub(ksub), V1(v1), max(max), min(min), Vector(AsVector) { }

    KN<T> Psub;
    KN<int> Ksub;
    KN<K> V1;
    K max;
    K min;
    bool Vector;
};
/**
 * @brief Class used to create a block of data which the server will send.
 * A basic ffPacket will look like this :
 * {
 *      "Version": "FreeFem++ Header 0.1",
 *      "Size": ** Compressed data size **
 * }
 * ** Header is padded using space to a size of 66 bytes using spaces **
 * ** The next data is compressed using the CBOR algorithm **
 * {
 *      "Geometry": [
 *          {
 *           "Type": ** Type of the data (eg: Curve, Mesh2D, Mesh3D) **,
 *           "ElementCount": ** Number of element in the array **,
 *           "ElementSize": ** Size of one element written as string formatted as "type number_of_elements ..." **,
 *           "Vertices": [...],
 *           "Indices": [...]
 *          },
 *          {},
 *          ...
 *      ]
 * }
 */
struct ffPacket {
    ffPacket();
    ffPacket(json data);

    ffPacket(const ffPacket& copy);
    ffPacket& operator=(const ffPacket& copy);

    std::string Dump(int indent = -1);
    void Compress();
    void Clear();

    std::string GetHeader();

    /**
     * @brief Convert the templated type to JSON
     *
     * @param const T& - Data to convert
     */
    template <typename T>
    void Jsonify(const T& data, long int Id);

    template <typename T>
    void JsonifyArgs(std::string Name, const T Data);

    json m_JSON;
    std::vector<uint8_t> m_Data;
    json m_Header;
};

#endif // FF_PACKET_HPP