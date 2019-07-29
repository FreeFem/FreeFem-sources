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
// SUMMARY : main file for TCP client
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Quentin Tessier
// E-MAIL  : tessier.quentin.pro@gmail.com

#include <unistd.h>
#include <limits.h>
#include "ffClient.h"

/**
 * @brief Construct a new ffClient from a io_context
 */
ffClient::ffClient(asio::io_context& io_context)
    : m_Socket(io_context), m_Timer(io_context)
{}

/**
 * @brief Start the client connection
 *
 * @param tcp::resolver::results_type - List of potentials connections
 */
void ffClient::Start(tcp::resolver::results_type endpoints)
{
    m_Endpoints = endpoints;
    StartConnect(m_Endpoints.begin());
}

/**
 * @brief Stop the client connection
 */
void ffClient::Stop()
{
    stopped = true;
    m_Socket.close();
    m_Timer.cancel();
}

/**
 * @brief Try one of the connection in the list
 *
 * @param tcp::resolver::results_type::iterator - Tested connection
 */
void ffClient::StartConnect(tcp::resolver::results_type::iterator endpoint_iter)
{
    if (endpoint_iter != m_Endpoints.end()) {
        m_Socket.async_connect(endpoint_iter->endpoint(), [=](std::error_code const & err)
        {
            HandleConnect(endpoint_iter, err);
        });
    } else {
        Stop();
    }
}

/**
 * @brief Check if the connection was successful, if not go to the next connection in the list
 *
 * @param tcp::resolver::results_type::iterator - Tested connection
 * @param const std::error_code& - Containt the return value of tcp::socket::async_connect
 */
void ffClient::HandleConnect(tcp::resolver::results_type::iterator endpoint_iter, std::error_code const & err)
{
    if (stopped)
        return;
    if (!m_Socket.is_open())
        StartConnect(++endpoint_iter);
    else if (err) {
        m_Socket.close();
        StartConnect(++endpoint_iter);
    } else {
        std::cout << "Connected !\n";
        StartRead();
        StartWrite();
    }
}

/**
 * @brief If the client is connected, send a heartbeat
 */
void ffClient::StartWrite()
{
    if (stopped)
        return;
    asio::async_write(m_Socket, asio::buffer("\n", 1), [=](std::error_code const & err, size_t /*  */)
    {
        HandleWrite(err);
    });
}

/**
 * @brief Check if writing the heartbeat message was successful,
 * send a new one after 15 seconds
 *
 * @param const std::error_code& - Containt the return value of asio::async_write
 */
void ffClient::HandleWrite(std::error_code const & err)
{
    if (stopped)
        return;
    if (!err) {
        m_Timer.expires_after(std::chrono::seconds(15));
        m_Timer.async_wait([=](std::error_code const & err)
        {
            StartWrite();
        });
    } else {
        std::cout << "Error on heartbeat: " << err.message() << "\n";
        Stop();
    }
}

/**
 * @brief Read the header of the packet sent
 */
void ffClient::StartRead()
{
    asio::async_read(m_Socket, asio::dynamic_buffer(m_HeaderBuffer, 66),
            [=](std::error_code const & err, size_t n)
            {
                HandleRead(err, n);
            });
}

/**
 * @brief Read the rest of the message from the header data
 *
 * @param const std::error_code& - Containt the return value of asio::async_read
 * @param std::size_t - Number of element read
 */
void ffClient::HandleRead(const std::error_code& err, std::size_t n)
{
    using json = nlohmann::json;
    if (stopped)
        return;
    if (!err) {
        json j = json::parse(m_HeaderBuffer);
        size_t data_size = j["Size"];
        std::error_code read_err;
        asio::read(m_Socket, asio::dynamic_buffer(m_DataBuffer, data_size), read_err);
        if (read_err) {
            std::cout << "Error on reading: " << err.message() << "\n";
            Stop();
            return;
        } else {
            m_HeaderBuffer.clear();
            m_DataBuffer.clear();
        }
    } else {
        std::cout << "Error on receive: " << err.message() << "\n";
        Stop();
    }
}

int main(int argc, char* argv[])
{
  try
  {
    if (argc != 3)
    {
      std::cerr << "Usage: client <host> <port>\n";
      return 1;
    }

    asio::io_context io_context;
    tcp::resolver r(io_context);
    ffClient c(io_context);

    c.Start(r.resolve(argv[1], argv[2]));

    io_context.run();
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 0;
}