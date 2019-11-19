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
// SUMMARY : TCP Server methods declaration
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Quentin Tessier
// E-MAIL  : tessier.quentin.pro@gmail.com

#include "ffServer.hpp"

// ffServer

/**
 * @brief Construct a new ffServer on the default port with at least 1 thread in the pool
 *
 * @param int ThreadNumber  Number of thread used in the pool (default value: 1)
 */
ffServer::ffServer(int ThreadNumber)
    : ThreadCount(ThreadNumber), m_ThreadPool(), m_ioService(), m_Acceptor(m_ioService, tcp::endpoint(tcp::v4(), 12345)), m_Connections(), m_Packets()
{}

/**
 * @brief Construct a new ffServer on a custom port with at least 1 thread in the pool
 *
 * @param short int Port    Custom port
 * @param int ThreadNumber  Number of thread used in the pool (default value: 1)
 */
ffServer::ffServer(short int Port, int ThreadNumber)
    : ThreadCount(ThreadNumber), m_ThreadPool(), m_ioService(), m_Acceptor(m_ioService, tcp::endpoint(tcp::v4(), Port)), m_Connections(), m_Packets()
{}

/**
 * @brief Start the server, invoking asio::io_service::run on every threads in the pool
 */
void ffServer::Start()
{
    std::cout << "Starting server.\n";
    AcceptConnection();
    for (int i = 0; i < ThreadCount; i += 1) {
        m_ThreadPool.emplace_back([=] { m_ioService.run(); });
    }
}

/**
 * @brief Stop the server, invoking asio::io_service::stop and joining every threads in the pool
 */
void ffServer::Stop()
{
    std::cout << "Stopping server.\n";
    m_ioService.stop();
    for (int i = 0; i < ThreadCount; i += 1) {
        m_ThreadPool[i].join();
    }
}

/**
 * @brief Wait for a new connection, invoking ffServer::HandleConnection on success
 */
void ffServer::AcceptConnection()
{
    ConnectHandle::pConnectHandle new_connect = ConnectHandle::create(m_ioService);
    m_Acceptor.async_accept(new_connect->m_Socket, std::bind(&ffServer::HandleConnection, this, new_connect, std::placeholders::_1));
}

/**
 * @brief Add the new connection to the list m_Connections if the connection was successful,
 * Calling ffServer::AcceptConnection creating a infinite loop of connection
 */
void ffServer::HandleConnection(ConnectHandle::pConnectHandle new_connect, std::error_code const & err)
{
    if (!err) {
        std::cout << "New Connection\n";
        m_Connections.push_back(new_connect);
        std::cout << "Reading client heartbeat.\n";
        if (!m_Packets.empty())
            SendAllStoredPacket(new_connect);
        StartRead(new_connect);
    }
    AcceptConnection();
}

/**
 * @brief Check if the last action on a client didn't timeout, removing the client from the list
 */
void ffServer::CheckDeadline(ConnectHandle::pConnectHandle connect, std::error_code const & /* */)
{
    if (connect->m_Deadline.expiry() <= asio::steady_timer::clock_type::now()) {
        connect->m_Socket.close();
        m_Connections.remove(connect);
        connect->m_Deadline.expires_at(asio::steady_timer::time_point::max());
    }
    connect->m_Deadline.async_wait([=](std::error_code const & err)
    {
        CheckDeadline(connect, err);
    });
}

/**
 * @brief Read data coming on the socket, this operation will timeout after 30 seconds
 */
void ffServer::StartRead(ConnectHandle::pConnectHandle connect)
{
    connect->m_Deadline.expires_after(std::chrono::seconds(30));

    connect->m_Deadline.async_wait([=](std::error_code const & err)
    {
        CheckDeadline(connect, err);
    });
    asio::async_read_until(connect->m_Socket, asio::dynamic_buffer(connect->m_MessageBuffer), '\n',
        [=](std::error_code const & err, size_t bytes_transferred)
        {
            HandleRead(connect, err, bytes_transferred);
        });
}

/**
 * @brief Check if asoi::async_read failed, invoking ffServer::StartRead making loop
 */
void ffServer::HandleRead(ConnectHandle::pConnectHandle connect, std::error_code const & err, size_t /* */)
{
    if (!err) {
        std::cout << "Client alive.\n";
        connect->m_MessageBuffer.clear();
    } else {
        std::cout << "Client disconnected\n";
        m_Connections.remove(connect);
        return;
    }
    StartRead(connect);
}

/**
 * @brief Send the packet to all connected clients
 */
void ffServer::Send(ffPacket packet)
{
    packet.m_JSON["Plot"] = m_Packets.size();
    packet.Compress();

    asio::steady_timer timer(m_ioService);
    size_t s = packet.m_Data.size();
    size_t i = s / MAX_PACKET_SIZE;

    for (auto & ite : m_Connections) {

        for (int j = 0; j < i; ++j) {
            timer.expires_after(std::chrono::seconds(1));
            packet.m_Header["Size"] = MAX_PACKET_SIZE;
            packet.m_Header["IDs"] = {m_Packets.size(), j};
            packet.m_Header["MaxPacket"] = i + ((s % MAX_PACKET_SIZE != 0) ? 1 : 0);
            ite->m_Socket.write_some(asio::buffer(packet.GetHeader(), MAX_PACKET_SIZE));

            ite->m_Socket.write_some(asio::buffer(packet.m_Data.data() + (MAX_PACKET_SIZE * j), MAX_PACKET_SIZE));
            //asio::write(tmp->m_Socket, asio::buffer(packet.m_Data));
            timer.wait();
        }
        if (s % MAX_PACKET_SIZE != 0) {
            packet.m_Header["Size"] = s % MAX_PACKET_SIZE;
            packet.m_Header["IDs"] = {m_Packets.size(), i + 1};
            packet.m_Header["MaxPacket"] = i + ((s % MAX_PACKET_SIZE != 0) ? 1 : 0);
            ite->m_Socket.write_some(asio::buffer(packet.GetHeader(), MAX_PACKET_SIZE));

            ite->m_Socket.write_some(asio::buffer(packet.m_Data.data() + (MAX_PACKET_SIZE * i), s % MAX_PACKET_SIZE));
        }
    }
    m_Packets.push_back(packet);
}

/**
 * @brief Send all the stored packet
 */
void ffServer::SendAllStoredPacket(ConnectHandle::pConnectHandle connect)
{
    auto ite = m_Packets.begin();
    asio::steady_timer timer(m_ioService);
    std::error_code err;

    for (ffPacket &p : m_Packets) {
        timer.expires_after(std::chrono::seconds(1));
        size_t s = p.m_Data.size();
        size_t i = s / MAX_PACKET_SIZE;

        for (auto & ite : m_Connections) {

            for (int j = 0; j < i; ++j) {
                p.m_Header["Size"] = MAX_PACKET_SIZE;
                p.m_Header["IDs"] = {m_Packets.size(), j};
                p.m_Header["MaxPacket"] = i + ((s % MAX_PACKET_SIZE != 0) ? 1 : 0);
                ite->m_Socket.write_some(asio::buffer(p.GetHeader(), MAX_PACKET_SIZE));

                ite->m_Socket.write_some(asio::buffer(p.m_Data.data() + (MAX_PACKET_SIZE * j), MAX_PACKET_SIZE));
                //asio::write(tmp->m_Socket, asio::buffer(packet.m_Data));
            }
            if (s % MAX_PACKET_SIZE != 0) {
                p.m_Header["Size"] = s % MAX_PACKET_SIZE;
                p.m_Header["IDs"] = {m_Packets.size(), i + 1};
                p.m_Header["MaxPacket"] = i + ((s % MAX_PACKET_SIZE != 0) ? 1 : 0);
                ite->m_Socket.write_some(asio::buffer(p.GetHeader(), MAX_PACKET_SIZE));

                ite->m_Socket.write_some(asio::buffer(p.m_Data.data() + (MAX_PACKET_SIZE * i), s % MAX_PACKET_SIZE));
            }
        }
        timer.wait();
        ite++;
    }
}

// ConnectHandle

ConnectHandle::ConnectHandle(asio::io_service& io_service)
    : m_Socket(io_service), m_Deadline(io_service.get_executor())
{}