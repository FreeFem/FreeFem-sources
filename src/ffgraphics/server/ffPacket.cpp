#include "ffPacket.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "Mesh3dn.hpp"

ffPacket::ffPacket()
{
    m_Header["Size"] = 0;
    m_Header["Version"] = "FreeFem++ Header 0.1";
    m_Data = std::vector<uint8_t>();
    m_JSON = json();
}

ffPacket::ffPacket(json data)
{
    m_Header["Size"] = 0;
    m_Header["Version"] = "FreeFem++ Header 0.1";

    m_JSON = json(data);
    m_Data = std::vector<uint8_t>();
}

ffPacket::ffPacket(const ffPacket& copy)
    : m_JSON(copy.m_JSON), m_Data(copy.m_Data), m_Header(copy.m_Header)
{}

ffPacket& ffPacket::operator=(const ffPacket& copy) {
    m_JSON = json(copy.m_JSON);
    m_Data = std::vector<uint8_t>(copy.m_Data);
    m_Header = json(copy.m_Header);

    return *this;
}

std::string ffPacket::Dump(int indent)
{
    return m_JSON.dump(indent);
}

void ffPacket::Compress()
{
    std::vector<uint8_t> tmp = json::to_cbor(m_JSON);

    m_Data.insert(m_Data.end(), tmp.begin(), tmp.end());
    m_Header["Size"] = m_Data.size();
}

void ffPacket::Clear()
{
    m_Header["Size"].clear();
    m_Data.clear();
    m_JSON.clear();
}

std::string ffPacket::GetHeader()
{
    std::string ret = m_Header.dump();
    while (ret.length() < 66)
        ret += ' ';
    return (ret);
}

/**
 * Transform a number of KN_<double> into JSON.
 * Sample JSON file :
 * {
 *      "Curve": {
 *          "_comment": "Dimension describe the dimension of the curve and how it's stored in JSON"
 *          "Dimension":
 *          "_comment": "Value is a array of double, where all data are stored in a contagious array ordered by point"
 *          "Values": []
 *      }
 * }
 */
template<>
void ffPacket::Jsonify<std::vector<KN_<double>>>(const std::vector<KN_<double>>& data)
{
    // Store the number of dimension for the curve (Max : 4)
    m_JSON["Curve"]["Dimension"] = data.size();
    int element_nbr = data[0].size(); // Number of sampled point on the function

    // For every sampled point store the double in the array "Values"
    for (int i = 0; i < element_nbr; i += 1) {
        // Store each dimension, add the value at index to the array
        for (int j = 0; j < 4; j += 1) {
            // If in range of number in dimension add a double else add a filler value
            if (j < data.size()) {
                m_JSON["Curve"]["Values"] += data[j][i];
            } else {
                m_JSON["Curve"]["Values"] += 0.f;
            }
        }
    }
}

template<>
void ffPacket::Jsonify<Fem2D::Mesh>(const Fem2D::Mesh& data)
{
    // Extract all vertices off the mesh
    for (int i = 0; i < data.nv; i += 1) {
        const Fem2D::Mesh::Vertex& p = data(i);
        m_JSON["Mesh"]["Vertices"] += p.x;
        m_JSON["Mesh"]["Vertices"] += p.y;
    }
    // Extract all elements off the mesh
    for (int i = 0; i < data.nt; i += 1) {
        const Fem2D::Mesh::Element& e(data[i]);
        m_JSON["Mesh"]["Elements"] += data(e[0]);
        m_JSON["Mesh"]["Elements"] += data(e[1]);
        m_JSON["Mesh"]["Elements"] += data(e[2]);
        m_JSON["Mesh"]["Elements"] += e.lab;
    }
    // Extract all borders off the mesh
    for (int i = 0; i < data.nbBrdElmts(); i += 1) {
        const Fem2D::Mesh::BorderElement& b(data.be(i));
        m_JSON["Mesh"]["Borders"] += data(b[0]);
        m_JSON["Mesh"]["Borders"] += data(b[1]);
        m_JSON["Mesh"]["Borders"] += b.lab;
    }
}

template<>
void ffPacket::Jsonify<Fem2D::Mesh3>(const Fem2D::Mesh3& data)
{
    // Extract all vertices off the mesh
    for (int i = 0; i < data.nv; i += 1) {
        const Fem2D::TVertex<Fem2D::R3> p = data.vertices[i];
        m_JSON["Mesh3"]["Vertices"] += p.x;
        m_JSON["Mesh3"]["Vertices"] += p.y;
        m_JSON["Mesh3"]["Vertices"] += p.z;
        m_JSON["Mesh3"]["Vertices"] += p.lab;
    }
    // Extract all elements off the mesh
    for (int i = 0; i < data.nt; i += 1) {
        const auto& e(data.elements[i]);
        m_JSON["Mesh3"]["Elements"] += data(e[0]);
        m_JSON["Mesh3"]["Elements"] += data(e[1]);
        m_JSON["Mesh3"]["Elements"] += data(e[2]);
        m_JSON["Mesh3"]["Elements"] += data(e[3]);
        m_JSON["Mesh3"]["Elements"] += e.lab;
    }
    // Extract all borders off the mesh
    for (int i = 0; i < data.nbBrdElmts(); i += 1) {
        const auto& b(data.borderelements[i]);
        m_JSON["Mesh3"]["Borders"] += data(b[0]);
        m_JSON["Mesh3"]["Borders"] += data(b[1]);
        m_JSON["Mesh3"]["Borders"] += data(b[2]);
        m_JSON["Mesh3"]["Borders"] += b.lab;
    }
}

template<>
void ffPacket::Jsonify<Fem2D::MeshS>(const Fem2D::MeshS& data)
{
    // Extract all vertices off the mesh
    for (int i = 0; i < data.nv; i += 1) {
        const Fem2D::TVertex<Fem2D::R3> p = data.vertices[i];
        m_JSON["Mesh3"]["Vertices"] += p.x;
        m_JSON["Mesh3"]["Vertices"] += p.y;
        m_JSON["Mesh3"]["Vertices"] += p.z;
        m_JSON["Mesh3"]["Vertices"] += p.lab;
    }
    // Extract all elements off the mesh
    for (int i = 0; i < data.nt; i += 1) {
        const auto& e(data.elements[i]);
        m_JSON["Mesh3"]["Elements"] += data(e[0]);
        m_JSON["Mesh3"]["Elements"] += data(e[1]);
        m_JSON["Mesh3"]["Elements"] += data(e[2]);
        m_JSON["Mesh3"]["Elements"] += data(e[3]);
        m_JSON["Mesh3"]["Elements"] += e.lab;
    }
    // Extract all borders off the mesh
    for (int i = 0; i < data.nbBrdElmts(); i += 1) {
        const auto& b(data.borderelements[i]);
        m_JSON["Mesh3"]["Borders"] += data(b[0]);
        m_JSON["Mesh3"]["Borders"] += data(b[1]);
        m_JSON["Mesh3"]["Borders"] += b.lab;
    }
}

template<>
void ffPacket::Jsonify<ffFE<Fem2D::R2, Fem2D::R>>(const ffFE<Fem2D::R2, Fem2D::R>& data)
{
    m_JSON["FE2D"]["Types"] = { "R2", "R" };

    for (int i = 0; i < data.Psub.N(); i += 1) {
        m_JSON["FE2D"]["Psub"] += data.Psub[i].x;
        m_JSON["FE2D"]["Psub"] += data.Psub[i].y;
    }
    for (int i = 0; i < data.Ksub.N(); i += 1) {
        m_JSON["FE2D"]["Ksub"] += data.Ksub[i];
    }
    for (int i = 0; i < data.V1.N(); i += 1) {
        m_JSON["FE2D"]["V1"] += data.V1[i];
    }
}

template<>
void ffPacket::Jsonify<ffFE<Fem2D::R2, complex<double>>>(const ffFE<Fem2D::R2, complex<double>>& data)
{
    m_JSON["FE2D"]["Types"] = { "R2", "complex<double>" };

    for (int i = 0; i < data.Psub.N(); i += 1) {
        m_JSON["FE2D"]["Psub"] += data.Psub[i].x;
        m_JSON["FE2D"]["Psub"] += data.Psub[i].y;
    }
    for (int i = 0; i < data.Ksub.N(); i += 1) {
        m_JSON["FE2D"]["Ksub"] += data.Ksub[i];
    }
    for (int i = 0; i < data.V1.N(); i += 1) {
        m_JSON["FE2D"]["V1"] += data.V1[i].real();
        m_JSON["FE2D"]["V1"] += data.V1[i].imag();
    }
}

template<>
void ffPacket::Jsonify<ffFE<Fem2D::R3, Fem2D::R>>(const ffFE<Fem2D::R3, Fem2D::R>& data)
{
    m_JSON["FE2D"]["Types"] = { "R3", "R" };

    for (int i = 0; i < data.Psub.N(); i += 1) {
        m_JSON["FE2D"]["Psub"] += data.Psub[i].x;
        m_JSON["FE2D"]["Psub"] += data.Psub[i].y;
        m_JSON["FE2D"]["Psub"] += data.Psub[i].z;
    }
    for (int i = 0; i < data.Ksub.N(); i += 1) {
        m_JSON["FE2D"]["Ksub"] += data.Ksub[i];
    }
    for (int i = 0; i < data.V1.N(); i += 1) {
        m_JSON["FE2D"]["V1"] += data.V1[i];
    }
}

template<>
void ffPacket::Jsonify<ffFE<Fem2D::R3, complex<double>>>(const ffFE<Fem2D::R3, complex<double>>& data)
{
    m_JSON["FE3D"]["Types"] = { "R3", "complex<double>" };

    for (int i = 0; i < data.Psub.N(); i += 1) {
        m_JSON["FE3D"]["Psub"] += data.Psub[i].x;
        m_JSON["FE3D"]["Psub"] += data.Psub[i].y;
        m_JSON["FE3D"]["Psub"] += data.Psub[i].z;
    }
    for (int i = 0; i < data.Ksub.N(); i += 1) {
        m_JSON["FE3D"]["Ksub"] += data.Ksub[i];
    }
    for (int i = 0; i < data.V1.N(); i += 1) {
        m_JSON["FE3D"]["V1"] += data.V1[i].real();
        m_JSON["FE3D"]["V1"] += data.V1[i].imag();
    }
}