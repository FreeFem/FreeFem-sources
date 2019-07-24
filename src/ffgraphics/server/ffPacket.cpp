#include "ffPacket.hpp"
#include "RNM.hpp"

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
void ffPacket::jsonify<std::vector<KN_<double>>>(std::vector<KN_<double>> data)
{
    m_JSON["Curve"]["Dimension"] = data.size();
    int element_nbr = data[0].size();

    for (int i = 0; i < element_nbr; i += 1) {
        for (int j = 0; j < data.size(); j += 1) {
            m_JSON["Curve"]["Values"] += data[j][i];
        }
    }
}