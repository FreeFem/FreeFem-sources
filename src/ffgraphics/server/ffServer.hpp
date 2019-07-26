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

class Connection : public std::enable_shared_from_this<Connection> {
  public:
    typedef std::shared_ptr<Connection> pointer;

    static pointer create(asio::io_service& io_service) {
        return pointer(new Connection(io_service));
    }

    tcp::socket& socket() { return m_Socket; }
    std::string& message() { return m_Message; }

  private:
    Connection(asio::io_service& io_service) : m_Socket(io_service) {}

    tcp::socket m_Socket;
    std::string m_Message;
};

class ffServer {
  public:
    ffServer(int thread_number = 1)
        : m_ioService(),
          m_Acceptor(m_ioService, tcp::endpoint(tcp::v4(), 12345)),
          nThreadCount(thread_number) {}

    ~ffServer() {
        m_ioService.stop();
        for (int i = 0; i < nThreadCount; i += 1) {
            m_ThreadPool[i].join();
        }
    }

    void start() {
        std::cout << "Starting server\n";
        for (int i = 0; i < nThreadCount; i += 1) {
            m_ThreadPool.emplace_back([=] { start_accept(); m_ioService.run(); });
        }

    }

    void stop() {
        m_ioService.stop();
        for (int i = 0; i < nThreadCount; i += 1) {
            m_ThreadPool[i].join();
        }
    }

    void send(ffPacket packet) {
        auto ite = m_Connection.begin();
        m_Packets.push_back(packet);

        std::cout << "Broadcasting json.\n";
        while (ite != m_Connection.end()) {
            auto tmp = *ite;
            asio::async_write(tmp->socket(), asio::buffer(packet.m_Header.dump()),
                              [this, ite](std::error_code const& err,
                                          size_t bytes_transferred) {
                                  if (err) {
                                      m_Connection.erase(ite);
                                  } else {
                                      std::cout << "Sent header\n";
                                  }
                              });
            asio::async_write(tmp->socket(), asio::buffer(packet.m_Data),
                              [this, ite](std::error_code const& err,
                                          size_t bytes_transferred) {
                                  if (err) {
                                      m_Connection.erase(ite);
                                  } else {
                                      std::cout << "Sent Data\n";
                                  }
                              });
            ite++;
        }
    }

  private:
    void handle_accept(Connection::pointer new_connection,
                       const std::error_code& err) {
        if (!err) {
            m_Connection.push_back(new_connection);
            start_read(new_connection);
        }
        start_accept();
    }

    void start_read(Connection::pointer new_connection) {
        asio::steady_timer timer(m_ioService);
        bool timer_result = false;
        bool read_result = false;

        timer.expires_after(std::chrono::seconds(10));
        timer.async_wait([&timer_result](std::error_code const & err) {
            timer_result = true;
        });
        asio::async_read(new_connection->socket(), asio::buffer(new_connection->message()),
                        [&read_result](std::error_code const & err, size_t /* */) {
                            read_result = true;
                        });
        while (!read_result && !timer_result) {
            if (read_result)
                timer.cancel();
            else if (timer_result) {
                std::cout << "Removing a client.\n";
                m_Connection.remove(new_connection);
            }
        }
        start_read(new_connection);
    }

    void start_accept() {
        std::cout << "start_accept\n";
        Connection::pointer new_connection = Connection::create(m_ioService);

        m_Acceptor.async_accept(new_connection->socket(),
                                std::bind(&ffServer::handle_accept, this,
                                          new_connection,
                                          std::placeholders::_1));
    }

    int nThreadCount;
    std::vector<std::thread> m_ThreadPool;
    asio::io_service m_ioService;
    tcp::acceptor m_Acceptor;
    std::list<Connection::pointer> m_Connection;
    std::list<ffPacket> m_Packets;
};

extern ffServer *graphicServer;

#endif // FFSERVER_HPP
