/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : TCP Client class
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Quentin Tessier
// E-MAIL  : tessier.quentin.pro@gmail.com

#include <asio/buffer.hpp>
#include <asio/io_context.hpp>
#include <asio/ip/tcp.hpp>
#include <asio/read_until.hpp>
#include <asio/steady_timer.hpp>
#include <asio/write.hpp>
#include <nlohmann/json.hpp>
#include <functional>
#include <iostream>
#include <string>
#include <deque>
#include <asio/read.hpp>
#include <utility>

using asio::steady_timer;
using asio::ip::tcp;
using std::placeholders::_1;
using std::placeholders::_2;

class ffClient {
    private:
        tcp::socket m_Socket;
        std::string m_HeaderBuffer;
        std::vector<uint8_t> m_DataBuffer;
        asio::steady_timer m_Timer;
        tcp::resolver::results_type m_Endpoints;

    public:
        bool stopped = false;

        ffClient(asio::io_context& io_context);

        void Start(tcp::resolver::results_type endpoints);
        void Stop();

        void StartConnect(tcp::resolver::results_type::iterator endpoint_iter);
        void HandleConnect(tcp::resolver::results_type::iterator endpoint_iter, std::error_code const & err);

        void StartWrite();
        void HandleWrite(const std::error_code& err);

        void StartRead();
        void HandleRead(const std::error_code& error, std::size_t n);
};