#ifndef  MODE_OPEN_HPP
#define  MODE_OPEN_HPP
#ifdef __WIN32__
#define  MODE_READ_BINARY "rb"
#define  MODE_WRITE_BINARY "wb"
#else
#define  MODE_READ_BINARY "r"
#define  MODE_WRITE_BINARY "w"
#endif
#endif
