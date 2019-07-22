#ifndef FF_PACKET_HPP
#define FF_PACKET_HPP

#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

struct ffPacket {
    std::string m_Header;
    std::vector<uint8_t> m_Data;

    void display_content() {
        std::cout << m_Header << "\n";
        for (size_t i = 0; i < m_Data.size(); ++i) {
            std::cout << m_Data[i];
        }
        std::cout << "\n";
    }
};

class ffSerializer {
    using json = nlohmann::json; // Quality of life using

  public:
    /* Build the basic header data */
    ffSerializer() {
        m_Header["Size"] = 0;
        m_Header["Version"] = "FreeFem++ Header 0.1";
        m_Header["Padding"] = "";
    }

    ffPacket CreatePacket(json data) {
        // Compress the JSON using the cbor algorithm
        std::vector<uint8_t> compressed = json::to_cbor(data);
        size_t size = compressed.size();
        std::string header;
        ffPacket packet;

        // Adding the real size of the message to the Header
        m_Header["Size"] = size;
        // Computing the padding necessary to have a Header of 66 bytes
        header = m_Header.dump();
        if (header.length() < 66) {
            std::string padding = "";
            for (size_t i = 0; i < 66 - header.length(); ++i)
                padding.push_back('a');
            m_Header["Padding"] = padding;
        }
        // Create the packet from it's header and packet
        packet.m_Data = compressed;
        packet.m_Header = m_Header.dump();
        packet.m_Header += "\n";
        return packet;
    }

    void AddData(json::value_type& vType, double data) { vType += data; }

  private:
    json m_Header;
};

#endif // FF_PACKET_HPP