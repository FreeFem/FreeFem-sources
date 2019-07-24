#ifndef FF_PACKET_HPP
#define FF_PACKET_HPP

#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;

class ffPacket {
    public:

        inline std::string& GetPacketHeader() { return m_Header; }

        inline std::vector<uint8_t>& GetPacketData() { return m_Data; }

        template <typename T>
        void jsonify(T data);

        inline std::string dump(int indent = -1) { return m_JSON.dump(indent); }

    private:

        void GetPacketHeader(int size)
        {
            json j;
            std::string tmp;

            j["Size"] = size;
            j["Version"] = "FreeFem++ Header 0.1";
            j["Padding"] = "";
            tmp = j.dump();
            if (tmp.length() < 66) {
                std::string padding = "";
                for (size_t i = 0; i < 66 - tmp.length(); ++i)
                    padding.push_back('a');
                j["Padding"] = padding;
            }
            m_Header = j.dump();
        }

        std::string m_Header;
        std::vector<uint8_t> m_Data;
        json m_JSON;
};

#endif // FF_PACKET_HPP