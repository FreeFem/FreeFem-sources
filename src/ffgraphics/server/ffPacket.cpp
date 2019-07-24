#include "ffPacket.hpp"
#include "RNM.hpp"


template<>
void ffPacket::jsonify<std::vector<KN_<double>>>(std::vector<KN_<double>> data)
{
    std::string names[4] = {"x", "y", "z", "w"};

    for (int i = 0; i < data.size(); i += 1) {
        double *tmp = data[i].operator double *();
        for (int j = 0; j < data[i].size(); j += 1)
            m_JSON["Curve"][names[i]] += tmp[j];
    }
}