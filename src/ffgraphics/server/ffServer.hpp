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

using asio::ip::tcp; // Quality of life using, get every asio's tcp class

class Connection
    : public std::enable_shared_from_this<Connection>
{
    public:
        typedef std::shared_ptr<Connection> pointer;

        static pointer create(asio::io_service& io_service)
        {
            return pointer(new Connection(io_service));
        }

        tcp::socket& socket() { return m_Socket; }

        void start()
        {
            m_Message = "Message\n\0";

            asio::async_write(m_Socket, asio::buffer(m_Message),
                    [=](std::error_code const & err, size_t bytes_transferred)
                    {
                        std::cout << "Send message of " << bytes_transferred << " bytes.\n";
                    });
        }

    private:
        Connection(asio::io_service& io_service)
            : m_Socket(io_service)
        { }

        tcp::socket m_Socket;
        std::string m_Message;
};

class ffServer {
    public:
        ffServer(asio::io_service& io_service)
            : m_ioService(io_service), m_Acceptor(io_service, tcp::endpoint(tcp::v4(), 12345))
        {
            start_accept();
        }

        void send(std::string message)
        {
            auto ite = m_Connection.begin();

            std::cout << "Broadcasting json.\n";
            while (ite != m_Connection.end()) {
                auto tmp = *ite;
                asio::async_write(tmp->socket(), asio::buffer(message),
                    [this, ite](std::error_code const & err, size_t bytes_transferred)
                    {
                        if (err) {
                            m_Connection.erase(ite);
                        } else {
                            std::cout << "Send message of " << bytes_transferred << " bytes.\n";
                        }
                    });
                ite++;
            }
        }

    private:
        void handle_accept(Connection::pointer new_connection, const std::error_code & err)
        {
            if (!err) {
                m_Connection.push_back(new_connection);
                new_connection->start();
            }
            start_accept();
        }

        void start_accept()
        {
            Connection::pointer new_connection = Connection::create(m_ioService);

            m_Acceptor.async_accept(new_connection->socket(),
                std::bind(&ffServer::handle_accept, this, new_connection, std::placeholders::_1));
        }

        asio::io_service& m_ioService;
        tcp::acceptor m_Acceptor;
        std::list<Connection::pointer> m_Connection;
};

#endif // FFSERVER_HPP
