#include "ffPacket.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "Mesh3dn.hpp"

static std::vector<float> ComputeColor(float compressed_color)
{
    float H = fabs(compressed_color) * 359.f;
    float S = 1.f;
    float V = 1.f;

    float C = V * S;
    float X = C * (1 - fabs(fmod(H / 60.f, 2.f) - 1));
    float m = V - C;

    float R, G, B = 0.f;
    if (H >= 0 && H < 60) {
        R = C; G = X; B = 0.f;
    } else if (H >= 60 && H < 120) {
        R = X; G = C; B = 0.f;
    } else if (H >= 120 && H < 180) {
        R = 0.f; G = C; B = X;
    } else if (H >= 180 && H < 240) {
        R = 0.f; G = X; B = C;
    } else if (H >= 240 && H < 300) {
        R = X; G = 0.f; B = C;
    } else if (H >= 300 && H < 360) {
        R = C; G = 0.f; B = X;
    }
    return {fabs((R + m)), fabs((G + m)), fabs((B + m)), 1.f};
}

static bool FindLabel(std::vector<long int> labels, long int nLabel)
{
    for (long int lab : labels) {
        if (lab == nLabel)
            return false;
    }
    return true;
}

static std::array<float, 4> LabelToColor(long lab)
{
    float H = 359.f / ((lab == 0) ? 1 : lab);
    float S = 1.f;
    float V = 1.f;

    float C = V * S;
    float X = C * (1 - fabs(fmod(H / 60.f, 2.f) - 1));
    float m = V - C;

    float R, G, B = 0.f;
    if (H >= 0 && H < 60) {
        R = C; G = X; B = 0.f;
    } else if (H >= 60 && H < 120) {
        R = X; G = C; B = 0.f;
    } else if (H >= 120 && H < 180) {
        R = 0.f; G = C; B = X;
    } else if (H >= 180 && H < 240) {
        R = 0.f; G = X; B = C;
    } else if (H >= 240 && H < 300) {
        R = X; G = 0.f; B = C;
    } else if (H >= 300 && H < 360) {
        R = C; G = 0.f; B = X;
    }
    return {fabs((R + m)), fabs((G + m)), fabs((B + m)), 1.f};
}

static int FindLabelIndex(std::vector<long int> ulabels, long int value)
{
    for (int i = 0; i < ulabels.size(); ++i) {
        if (ulabels[i] == value)
            return i;
    }
    return 0;
}

static std::array<float, 4> ComputeLabel(int ColorIndex, int ColorCount)
{
    long int H = 359 / ColorCount * ColorIndex;

    return LabelToColor(H);
}

/**
 * @brief Construct a new ffPacket
 */
ffPacket::ffPacket()
{
    m_Header["Size"] = 0;
    m_Header["Version"] = "FreeFem++ Header 0.1";
    m_Data = std::vector<uint8_t>();
    m_JSON = json();
    m_JSON["Geometry"] = json::array();
}

/**
 * @brief Construct a new ffPacket from JSON
 *
 * @param[in] json data - JSON data
 */
ffPacket::ffPacket(json data)
{
    m_Header["Size"] = 0;
    m_Header["Version"] = "FreeFem++ Header 0.1";

    m_JSON = json(data);
    m_JSON["Geometry"] = json::array();
    m_Data = std::vector<uint8_t>();
}

/**
 * @brief Copy constructor
 *
 * @param[in] const ffPacket& - packet to copy
 */
ffPacket::ffPacket(const ffPacket& copy)
    : m_JSON(copy.m_JSON), m_Data(copy.m_Data), m_Header(copy.m_Header)
{}

/**
 * @brief Copy operator
 *
 * @param[in] const ffPacket& - packet to copy
 * @return ffPacket& - return a new packet
 */
ffPacket& ffPacket::operator=(const ffPacket& copy) {
    m_JSON = json(copy.m_JSON);
    m_Data = std::vector<uint8_t>(copy.m_Data);
    m_Header = json(copy.m_Header);

    return *this;
}

/**
 * @brief Make a string out of the content of the packet
 *
 * @param[in] int - number of spaces used to indent (def: -1)
 * @return std::string - string
 */
std::string ffPacket::Dump(int indent)
{
    return m_JSON.dump(indent);
}

/**
 * @brief Compress the content of the packet and updating the header with the new size
 */
void ffPacket::Compress()
{
    //std::cout << Dump() << "\n";
    std::vector<uint8_t> tmp = json::to_cbor(m_JSON);

    m_Data.insert(m_Data.end(), tmp.begin(), tmp.end());
    m_Header["Size"] = m_Data.size();
}

/**
 * @brief Empty the packet
 */
void ffPacket::Clear()
{
    m_Header["Size"].clear();
    m_Data.clear();
    m_JSON.clear();
    m_JSON["Geometry"] += json::object();
}

/**
 * @brief Returns the header
 *
 * @return std::string - the header with the right number of octet
 */
std::string ffPacket::GetHeader()
{
    std::string ret = m_Header.dump();
    while (ret.length() < 64)
        ret += ' ';
    return (ret);
}

/**
 * @brief Specification of ffPacket::Jsonify for a curve
 * The JSON format for a curve is :
 *
 *  {
 *      "Curve": {
 *          "Dimension": n,
 *          "Values": [x0, y0, z0, w0, x1, y1, z1, w1, ...]
 *      }
 *  }
 *
 * @param const std::vector<KN_<double>>& - a curve
 */
template<>
void ffPacket::Jsonify<std::vector<KN_<double>>>(const std::vector<KN_<double>>& data)
{
    m_JSON["Geometry"] += json::object();
    auto Obj_JSON = --(m_JSON["Geometry"].end());


    (*Obj_JSON)["Primitive"] = "Line";
    (*Obj_JSON)["GeometryType"] = "Surface";
    (*Obj_JSON)["Vertices"] = json::array();
    (*Obj_JSON)["Indices"] = json::array();
    (*Obj_JSON)["Labels"] = json::array();
    size_t element_nbr = data[0].size();
    size_t dimensions = (data.size() == 4) ? 3 : data.size();

    if (data.size() == 2)
        (*Obj_JSON)["Type"] = "Mesh2D";
    else
        (*Obj_JSON)["Type"] = "Mesh3D";

    for (size_t i = 0; i < element_nbr; ++i) {
        for (size_t j = 0; j < dimensions; ++j) {
            (*Obj_JSON)["Vertices"] += data[j][(long long int)i];
        }
        if (data.size() == 4) {
            (*Obj_JSON)["Labels"] += (long int)data[3][(long long int)i];
        } else {
            (*Obj_JSON)["Labels"] += 0;
        }
        if (i != element_nbr - 1) {
            (*Obj_JSON)["Indices"] += i;
            (*Obj_JSON)["Indices"] += i + 1;
        }
    }
}

/**
 * @brief Specification of ffPacket::Jsonify for a Mesh
 * The JSON format for a Mesh is :
 *
 *  {
 *      "Mesh": {
 *          "Vertices": [x0, y0, x1, y1 ...],
 *          "Elements": [x0, y0, z0, w0, x1, y1, z1, w1, ...],
 *          "Borders": [x0, y0, z0, x1, y1, z1, ...]
 *      }
 *  }
 *
 * @param const Fem2D::Mesh& - a Mesh
 */
template<>
void ffPacket::Jsonify<Fem2D::Mesh>(const Fem2D::Mesh& data)
{
    m_JSON["Geometry"] += json::object();
    auto Obj_JSON = --(m_JSON["Geometry"].end());

    (*Obj_JSON)["Type"] = "Mesh2D";
    (*Obj_JSON)["GeometryType"] = "Volume";
    (*Obj_JSON)["Primitive"] = "Triangle";
    (*Obj_JSON)["Vertices"] = json::array();
    (*Obj_JSON)["Indices"] = json::array();
    (*Obj_JSON)["Labels"] = json::array();

    std::cout << "Number of triangles : " << data.nt << "\n";
    for (int i = 0; i < data.nt; ++i) {
        const Fem2D::Mesh::Element& e(data[i]);
        (*Obj_JSON)["Indices"] += data(e[0]);
        (*Obj_JSON)["Indices"] += data(e[1]);
        (*Obj_JSON)["Indices"] += data(e[2]);

        (*Obj_JSON)["Labels"] += e.lab;
        (*Obj_JSON)["Labels"] += e.lab;
        (*Obj_JSON)["Labels"] += e.lab;
    }
    for (int i = 0; i < data.nv; ++i) {
        const Fem2D::Mesh::Vertex& p = data(i);
        (*Obj_JSON)["Vertices"] += p.x;
        (*Obj_JSON)["Vertices"] += p.y;
    }

    // m_JSON["Geometry"] += json::object();
    // Obj_JSON = --(m_JSON["Geometry"].end());

    // (*Obj_JSON)["Type"] = "Mesh2D";
    // (*Obj_JSON)["GeometryType"] = "Surface";
    // (*Obj_JSON)["Primitive"] = "Line";
    // (*Obj_JSON)["Vertices"] = json::array();
    // (*Obj_JSON)["Indices"] = json::array();

    // int index = 0;
    // for (int i = 0; i < data.nbBrdElmts(); ++i) {
    //     const Fem2D::Mesh::BorderElement& b(data.be(i));
    //     const Fem2D::Mesh::Vertex& p = data(data(b[0]));
    //     const Fem2D::Mesh::Vertex& p1 = data(data(b[1]));

    //     (*Obj_JSON)["Vertices"] += p.x;
    //     (*Obj_JSON)["Vertices"] += p.y;
    //     (*Obj_JSON)["Vertices"] += 0.0f;
    //     std::array<float, 4> Color = LabelToColor(b.lab);
    //     (*Obj_JSON)["Vertices"] += Color[0];
    //     (*Obj_JSON)["Vertices"] += Color[1];
    //     (*Obj_JSON)["Vertices"] += Color[2];
    //     (*Obj_JSON)["Vertices"] += Color[3];

    //     (*Obj_JSON)["Indices"] += index;
    //     index += 1;

    //     (*Obj_JSON)["Vertices"] += p1.x;
    //     (*Obj_JSON)["Vertices"] += p1.y;
    //     (*Obj_JSON)["Vertices"] += 0.0f;
    //     Color = LabelToColor(b.lab);
    //     (*Obj_JSON)["Vertices"] += Color[0];
    //     (*Obj_JSON)["Vertices"] += Color[1];
    //     (*Obj_JSON)["Vertices"] += Color[2];
    //     (*Obj_JSON)["Vertices"] += Color[3];

    //     (*Obj_JSON)["Indices"] += index;
    //     index += 1;
    // }
}

static void PushVertex(json& VerticesArray, std::array<std::array<float, 21>, 4> triangles)
{
    for (int i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < triangles[i].size(); ++j) {
            VerticesArray += triangles[i][j];
        }
    }
}

/**
 * @brief Specification of ffPacket::Jsonify for a Mesh3
 * The JSON format for a Mesh3 is :
 *
 *  {
 *      "Mesh3": {
 *          "Vertices": [x0, y0, z0, w0, x1, y1, z1, w1 ...],
 *          "Elements": [x0, y0, z0, w0, v0, x1, y1, z1, w1, v0, ...],
 *          "Borders": [x0, y0, z0, w0, x1, y1, z1, w1, ...]
 *      }
 *  }
 *
 * @param const Fem2D::Mesh3& - a Mesh3
 */
template<>
void ffPacket::Jsonify<Fem2D::Mesh3>(const Fem2D::Mesh3& data)
{
    m_JSON["Geometry"] += json::object();
    auto Obj_JSON = --(m_JSON["Geometry"].end());

    (*Obj_JSON)["Type"] = "Mesh3D";
    (*Obj_JSON)["GeometryType"] = "Volume";
    (*Obj_JSON)["Primitive"] = "Triangle";
    (*Obj_JSON)["Vertices"] = json::array();
    (*Obj_JSON)["Indices"] = json::array();
    (*Obj_JSON)["Labels"] = json::array();

    for (int i = 0; i < data.nv; ++i) {
        (*Obj_JSON)["Vertices"] += data.vertices[i].x;
        (*Obj_JSON)["Vertices"] += data.vertices[i].y;
        (*Obj_JSON)["Vertices"] += data.vertices[i].z;
    }
    for (int i = 0; i < data.nt; ++i) {
        const auto& tetra(data.elements[i]);
        // Triangle 1
        (*Obj_JSON)["Indices"] += data(tetra[0]);
        (*Obj_JSON)["Indices"] += data(tetra[1]);
        (*Obj_JSON)["Indices"] += data(tetra[2]);
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;

        // Triangle 2
        (*Obj_JSON)["Indices"] += data(tetra[0]);
        (*Obj_JSON)["Indices"] += data(tetra[2]);
        (*Obj_JSON)["Indices"] += data(tetra[3]);
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;

        // Triangle 3
        (*Obj_JSON)["Indices"] += data(tetra[0]);
        (*Obj_JSON)["Indices"] += data(tetra[3]);
        (*Obj_JSON)["Indices"] += data(tetra[1]);
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;

        // Triangle 4
        (*Obj_JSON)["Indices"] += data(tetra[1]);
        (*Obj_JSON)["Indices"] += data(tetra[2]);
        (*Obj_JSON)["Indices"] += data(tetra[3]);
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;
        (*Obj_JSON)["Labels"] += tetra.lab;
    }

    m_JSON["Geometry"] += json::object();
    Obj_JSON = --(m_JSON["Geometry"].end());

    (*Obj_JSON)["Type"] = "Mesh3D";
    (*Obj_JSON)["GeometryType"] = "Surface";
    (*Obj_JSON)["Primitive"] = "Triangle";
    (*Obj_JSON)["Vertices"] = json::array();
    (*Obj_JSON)["Indices"] = json::array();
    (*Obj_JSON)["Labels"] = json::array();

    int count = 0;
    for (int i = 0; i < data.nbe; i += 1) {
        const auto& b(data.borderelements[i]);
        (*Obj_JSON)["Vertices"] += b[0].x;
        (*Obj_JSON)["Vertices"] += b[0].y;
        (*Obj_JSON)["Vertices"] += b[0].z;

        (*Obj_JSON)["Vertices"] += b[1].x;
        (*Obj_JSON)["Vertices"] += b[1].y;
        (*Obj_JSON)["Vertices"] += b[1].z;

        (*Obj_JSON)["Vertices"] += b[2].x;
        (*Obj_JSON)["Vertices"] += b[2].y;
        (*Obj_JSON)["Vertices"] += b[2].z;

        (*Obj_JSON)["Indices"] += count;
        (*Obj_JSON)["Indices"] += count + 1;
        (*Obj_JSON)["Indices"] += count + 2;
        count += 3;
        (*Obj_JSON)["Labels"] += b.lab;
        (*Obj_JSON)["Labels"] += b.lab;
        (*Obj_JSON)["Labels"] += b.lab;
    }
}

/**
 * @brief Specification of ffPacket::Jsonify for a MeshS
 * The JSON format for a MeshS is :
 *
 *  {
 *      "MeshS": {
 *          "Vertices": [x0, y0, z0, w0, x1, y1, z1, w1 ...],
 *          "Elements": [x0, y0, z0, w0, v0, x1, y1, z1, w1, v0, ...],
 *          "Borders": [x0, y0, z0, x1, y1, z1, ...]
 *      }
 *  }
 *
 * @param const Fem2D::MeshS& - a MeshS
 */
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

/**
 * @brief Specification of ffPacket::Jsonify for a 2D FEspace using real
 * The JSON format for a 2D FEspace using real is :
 *
 *  {
 *      "FE2D": {
 *          "Type": { "R2", "R" },
 *          "Psub": [x0, y0, x1, y1, ...],
 *          "Ksub": [x0, x1, ...],
 *          "V1": [x0, x1, ...]
 *      }
 *  }
 *
 * @param const ffFE<Fem2D::R2, Fem2D::R>& - a 2D FEspace using real
 */
template<>
void ffPacket::Jsonify<ffFE<Fem2D::R2, Fem2D::R>>(const ffFE<Fem2D::R2, Fem2D::R>& data)
{
    m_JSON["Geometry"] += json::object();
    auto Obj_JSON = --(m_JSON["Geometry"].end());

    (*Obj_JSON)["Type"] = "Mesh2D";
    (*Obj_JSON)["GeometryType"] = "Surface";
    (*Obj_JSON)["Primitive"] = "Line";
    (*Obj_JSON)["Vertices"] = json::array();
    (*Obj_JSON)["Indices"] = json::array();
    (*Obj_JSON)["Labels"] = json::array();

    for (int i = 0; i < data.Psub.N(); i += 1) {
        (*Obj_JSON)["Vertices"] += data.Psub[i].x;
        (*Obj_JSON)["Vertices"] += data.Psub[i].y;
    }
    for (int i = 0; i < data.Ksub.N(); i += 1) {
        m_JSON["FE2D"]["Ksub"] += data.Ksub[i];
    }
    for (int i = 0; i < data.V1.N(); i += 1) {
        m_JSON["FE2D"]["V1"] += data.V1[i];
    }
}

/**
 * @brief Specification of ffPacket::Jsonify for a 2D FEspace using complex
 * The JSON format for a 2D FEspace using complex is :
 *
 *  {
 *      "FE2D": {
 *          "Type": { "R2", "complex<double>" },
 *          "Psub": [x0, y0, x1, y1, ...],
 *          "Ksub": [x0, x1, ...],
 *          "V1": [x0, x1, ...]
 *      }
 *  }
 *
 * @param const ffFE<Fem2D::R2, complex<double>>& - a 2D FEspace using complex
 */
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

/**
 * @brief Specification of ffPacket::Jsonify for a 3D FEspace using real
 * The JSON format for a 3D FEspace using real is :
 *
 *  {
 *      "FE2D": {
 *          "Type": { "R3", "R" },
 *          "Psub": [x0, y0, z0, x1, y1, z1, ...],
 *          "Ksub": [x0, x1, ...],
 *          "V1": [x0, x1, ...]
 *      }
 *  }
 *
 * @param const ffFE<Fem2D::R3, Fem2D::R>& - a 3D FEspace using real
 */
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

/**
 * @brief Specification of ffPacket::Jsonify for a 3D FEspace using complex
 * The JSON format for a 3D FEspace using complex is :
 *
 *  {
 *      "FE2D": {
 *          "Type": { "R3", "complex<double>" },
 *          "Psub": [x0, y0, z0, x1, y1, z1, ...],
 *          "Ksub": [x0, x1, ...],
 *          "V1": [x0, y0, x1, y1, ...]
 *      }
 *  }
 *
 * @param const ffFE<Fem2D::R3, complex<double>>& - a 3D FEspace using complex
 */
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

template<>
void ffPacket::JsonifyArgs<double>(std::string Name, const double Data)
{
    m_JSON[Name] += Data;
}

template<>
void ffPacket::JsonifyArgs<bool>(std::string Name, const bool Data)
{
    m_JSON[Name] += Data;
}

template<>
void ffPacket::JsonifyArgs<int>(std::string Name, const int Data)
{
    m_JSON[Name] += Data;
}

template<>
void ffPacket::JsonifyArgs<std::string>(std::string Name, const std::string Data)
{
    m_JSON[Name] += Data;
}

template<>
void ffPacket::JsonifyArgs<long int>(std::string Name, const long int Data)
{
    m_JSON[Name] += Data;
}

template<>
void ffPacket::JsonifyArgs<KN_<double>>(std::string Name, const KN_<double> Data)
{
    if (Data.size() == 0) {
        m_JSON[Name] += {0.f};
        return;
    }
    std::vector<double> Values(Data.size());

    for (size_t i = 0; i < Data.size(); ++i) {
        Values.push_back(Data[(long long int)i]);
    }
    m_JSON[Name] += Values;
}

template<>
void ffPacket::JsonifyArgs<KN<double>>(std::string Name, const KN<double> Data)
{
    if (Data.size() == 0) {
        m_JSON[Name] += {0.f};
        return;
    }
    std::vector<double> Values(Data.size());

    for (size_t i = 0; i < Data.size(); ++i) {
        Values.push_back(Data[(long long int)i]);
    }
    m_JSON[Name] += Values;
}