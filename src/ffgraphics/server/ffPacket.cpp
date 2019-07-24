#include "ffPacket.hpp"
#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"

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
void ffPacket::jsonify<std::vector<KN_<double>>>(const std::vector<KN_<double>>& data)
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
void ffPacket::jsonify<Fem2D::Mesh>(const Fem2D::Mesh& data)
{
    // Extract all vertices off the mesh
    for (int i = 0; i < data.nv; i += 1) {
        Fem2D::Mesh::Vertex& p = data(i);
        m_JSON["Mesh"]["Vertices"] += p.x;
        m_JSON["Mesh"]["Vertices"] += p.y;
    }
    // Extract all elements off the mesh
    for (int i = 0; i < data.nt; i += 1) {
        Fem2D::Mesh::Element& e(data[i]);
        m_JSON["Mesh"]["Elements"] += data(e[0]);
        m_JSON["Mesh"]["Elements"] += data(e[1]);
        m_JSON["Mesh"]["Elements"] += data(e[2]);
        m_JSON["Mesh"]["Elements"] += e.lab;
    }
    // Extract all borders off the mesh
    for (int i = 0; i < data.nbBrdElmts(); i += 1) {
        Fem2D::Mesh::BorderElement& b(data.be(i));
        m_JSON["Mesh"]["Borders"] += data(b[0]);
        m_JSON["Mesh"]["Borders"] += data(b[1]);
        m_JSON["Mesh"]["Borders"] += b.lab;
    }
}