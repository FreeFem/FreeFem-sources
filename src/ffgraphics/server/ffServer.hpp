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
// SUMMARY : TCP Server class
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Quentin Tessier
// E-MAIL  : tessier.quentin.pro@gmail.com


#ifndef FFSERVER_HPP
#define FFSERVER_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <deque>
#include <list>
#include <functional>
#include <asio.hpp>
#include "ffPacket.hpp"

using asio::ip::tcp; // Quality of life using, get every asio's tcp class

class ConnectHandle : public std::enable_shared_from_this<ConnectHandle>
{
    public:
        typedef std::shared_ptr<ConnectHandle> pConnectHandle;

        static pConnectHandle create(asio::io_service& io_service)
        { return pConnectHandle(new ConnectHandle(io_service)); }

        ConnectHandle(asio::io_service& io_service);

        tcp::socket m_Socket;
        asio::steady_timer m_Deadline;
        std::string m_MessageBuffer;
};

class ffServer {
    private:
        // ThreadPool on which the server will run
        int ThreadCount;
        std::vector<std::thread> m_ThreadPool;

        asio::io_service m_ioService;
        tcp::acceptor m_Acceptor;

        std::list<ConnectHandle::pConnectHandle> m_Connections;
        std::list<ffPacket> m_Packets;

    public:
        ffServer(int ThreadNumber = 1);
        ffServer(short int Port, int ThreadNumber = 1);

        // Deleting copy constructors (Using them will throw a error at compilation)
        ffServer(const ffServer& tocopy) = delete;
        ffServer& operator=(const ffServer& tocopy) = delete;

        void Start();
        void Stop();

        void AcceptConnection();
        void HandleConnection(ConnectHandle::pConnectHandle new_connect, std::error_code const & err);

        void CheckDeadline(ConnectHandle::pConnectHandle connect, std::error_code const & err);

        void StartRead(ConnectHandle::pConnectHandle connect);
        void HandleRead(ConnectHandle::pConnectHandle connect, std::error_code const & err, size_t bytes_transferred);

        void Send(ffPacket packet);
        void SendAllStoredPacket(ConnectHandle::pConnectHandle connect);
};

extern ffServer *graphicServer;

#endif // FFSERVER_HPP
