/****************************************************************************/
/* This file is part of FreeFem++.                                          */
/*                                                                          */
/* FreeFem++ is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFem++ is distributed in the hope that it will be useful,             */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        */
/****************************************************************************/
// SUMMARY : Tools to read ppm file
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

/* Usage in FreeFEM .edp file
 * see: examples/plugin/ppm2rnm.edp
 */
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;

#include "RNM.hpp"
#include <cmath>
typedef KNM<double> *pRnm;
typedef KN<double> *pRn;
typedef string *pstring;
#include "ppmimg.h"

#define DISP_INFO "PPM2RMN:"
#define DISP_ERROR "PPM2RNM - ERROR:"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Load a PPM image
 * \param imgname PPM image file name
 * \param type PPM image type
 * \param quiet Verbosity
 * \return Success: loaded PPM image, Error: 0
 */
PPMimage *loadPPM (const char *imgname, ubyte *type, ubyte quiet) {
	pPPMimage result;
	FILE *fp;
	int typimg, ret, s, maxval, bitsize;
	const char *ptr;
	char buff[1024], data[256];

	/* search for image */
	fprintf(stdout, "%s Loading image: %s\n", DISP_INFO, imgname);
	ptr = strstr(imgname, ".ppm");
	strcpy(data, imgname);
	if (!ptr) {
		ptr = strstr(imgname, ".pgm");
		if (!ptr) {strcat(data, ".ppm");}
	}

	fp = fopen(data, "rb");
	if (!fp) {
		fprintf(stderr, "%s UNABLE TO OPEN FILE %s.\n", DISP_ERROR, data);
		return 0;
	}

	if (!quiet)
		fprintf(stdout, "%s Opening %s\n", DISP_INFO, data);

	if (!fgets(buff, sizeof(buff), fp)) {
		fprintf(stderr, "%s INVALID HEADER.\n", DISP_ERROR);
		return 0;
	}

	/* check header file */
	if (buff[0] != 'P') {
		fprintf(stderr, "%s INVALID IMAGE FORMAT (MUST BE 'PX').\n", DISP_ERROR);
		return 0;
	}

	switch (buff[1]) {
		case '2': typimg = P2;
			break;
		case '3': typimg = P3;
			break;
		case '4': typimg = P4;
			break;
		case '5': typimg = P5;
			break;
		case '6': typimg = P6;
			break;
		default:
			fprintf(stderr, "%s INVALID IMAGE FORMAT (MUST BE 'PX').\n", DISP_ERROR);
			fclose(fp);
			return 0;
	}

	/* allocate memory to store image */
	result = (PPMimage *)malloc(sizeof(PPMimage));
	assert(result);

	do {
		ret = fscanf(fp, "%s", buff);
		if (ret == EOF) {break;}

		/* check and strip comments */
		if (buff[0] == '#') {
			char c;
			do {
				c = getc(fp);
			} while (c != '\n');
		} else {break;}
	} while (1);

	/* read columns + lines */
	ret = sscanf(buff, "%d", &s);
	result->sizeX = (short)s;
	ret += fscanf(fp, "%d", &s);
	result->sizeY = (short)s;
	if (ret != 2) {
		fprintf(stderr, "%s ERROR LOADING IMAGE.\n", DISP_ERROR);
		free(result);
		fclose(fp);
		return 0;
	}

	if (fscanf(fp, "%d", &maxval) != 1) {
		fprintf(stderr, "%s INVALID IMAGE SIZE.\n", DISP_ERROR);
		free(result);
		fclose(fp);
		return 0;
	}

	/* strip line */
	while (fgetc(fp) != '\n') {;}

	/* size based on type */
	if (typimg == P2 || typimg == P5 || typimg == P4) {
		bitsize = result->sizeX * result->sizeY;
	} else {
		bitsize = 3 * result->sizeX * result->sizeY;
	}

	if (!quiet) {
		fprintf(stdout, "%s Image size: %dx%d - %d bytes\n",
		        DISP_INFO, result->sizeX, result->sizeY, bitsize);
	}

	result->data = (ubyte *)malloc(1 + bitsize * sizeof(ubyte));
	assert(result->data);

	/* read data file */
	switch (typimg) {
	case P2:/* ascii file (grey)  */
	case P3:/* ascii file (color) */
		int i;
		for (i = 0; i < bitsize; i++) {
			int r;
			int rr = fscanf(fp, "%d", &r);
			if (rr == EOF) {
				fprintf(stderr, "%s ERROR LOADING IMAGE.\n", DISP_ERROR);
				free(result->data);
				free(result);
				fclose(fp);
				return 0;
			}
			result->data[i] = (ubyte)r;
		}
		break;

	case P5:/* binary file (grey) */
	case P6:/* binary file (color) */
		ret = fread(result->data, sizeof(ubyte), bitsize, fp);
		if (ret != bitsize) {
			fprintf(stderr, "%s ERROR LOADING IMAGE.\n", DISP_ERROR);
			free(result->data);
			free(result);
			fclose(fp);
			return 0;
		}
		break;
	}

	fclose(fp);

	if (*type == DEFAULT) {
		if (typimg == P2 || typimg == P5) {
			*type = GREY;
		} else {
			*type = COLOR;
		}
	}
	/* convert to grey levels */
	else if (*type == GREY && (typimg == P3 || typimg == P6)) {
		int i, k;
		fprintf(stdout, "%s Converting to grey levels\n", DISP_INFO);

		for (i = 0, k = 0; i < bitsize; i += 3, k++) {
			int r = (int)result->data[i];
			int g = (int)result->data[i + 1];
			int b = (int)result->data[i + 2];
			result->data[k] = (ubyte)(0.3 * r + 0.59 * g + 0.11 * b);
		}

		result->data = (ubyte *)realloc(result->data, sizeof(ubyte) * bitsize / 3 + 1);
	}

	return result;
}

/*!
 * \brief Save a PPM image
 * \param imgname PPM image file name
 * \param img PPM image
 * \param typimg PPM image type
 * \return Success: 1, Error: 0
 */
int savePPM (const char *imgname, pPPMimage img, int typimg) {
	FILE *out;
	int bitsize;

	/* open file */
	out = fopen(imgname, "w");
	if (!out) {
		fprintf(stderr, "%s UNABLE TO OPEN FILE %s.\n", DISP_ERROR, imgname);
		return 0;
	}

	/* write out image file */
	bitsize = img->sizeX * img->sizeY;

	switch (typimg) {
		case P2:
			int i, c;

			fprintf(out, "P2\n");
			fprintf(out, "# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
			fprintf(out, "%d %d\n", img->sizeX, img->sizeY);
			fprintf(out, "255\n");
			c = 0;

			for (i = 0; i < img->sizeX * img->sizeY; i++) {
				fprintf(out, "%3d ", (int)img->data[i]);
				if (++c == 17) {
					c = 0;
					fprintf(out, "\n");
				}
			}

			fprintf(out, "\n");
			break;
		case P5:
			fprintf(out, "P5\n");
			fprintf(out, "# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
			fprintf(out, "%d %d\n", img->sizeX, img->sizeY);
			fprintf(out, "255\n");
			fwrite(img->data, sizeof(ubyte), bitsize, out);
			break;
		case P6:
			fprintf(out, "P6\n");
			fprintf(out, "# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
			fprintf(out, "%d %d\n", img->sizeX, img->sizeY);
			fprintf(out, "255\n");
			fwrite(img->data, sizeof(ubyte), 3 * bitsize, out);
			break;
	}

	fclose(out);

	return 1;
}

/*! \brief Compute difference between two images
 * \param bits First PPM image
 * \param img Second PPM image
 * \param itype PPM image type
 * \return Success: Difference PPM image, Error: 0
 */
pPPMimage diffImg (pPPMimage bits, pPPMimage img, ubyte itype) {
	pPPMimage dif;
	double psnr;
	int i, bitsize, dmax;

	fprintf(stdout, "%s Image difference\n", DISP_INFO);
	bitsize = (int)bits->sizeX * bits->sizeY;
	if (itype == COLOR) {bitsize *= 3;}

	dif = (PPMimage *)malloc(sizeof(PPMimage));
	if (!dif) {
		fprintf(stderr, "%s Sorry, not enough memory. Bye.\n", DISP_ERROR);
		return 0;
	}

	dif->sizeX = bits->sizeX;
	dif->sizeY = bits->sizeY;
	dif->data = (ubyte *)malloc(bitsize * sizeof(ubyte));
	if (!dif->data) {
		fprintf(stderr, "%s Sorry, not enough memory. Bye.\n", DISP_ERROR);
		free(dif);
		return 0;
	}

	dmax = 0;
	psnr = 0.0f;

	for (i = 0; i < bitsize; i++) {
		double dd = abs((int)(bits->data[i] - img->data[i]));
		dmax = max(dmax, dd);
		psnr += (double)dd * dd;
		dif->data[i] = (ubyte)(255 - dd);
	}

	if (psnr == 0.0f) {fprintf(stderr, "%s PSNR problem!", DISP_ERROR);} else {
		psnr = 65025.0f / psnr;
		psnr = 10.0 * log10(bitsize * psnr);
	}

	fprintf(stdout, "%s PSNR = %.2f - dmax = %d\n", DISP_INFO, psnr, dmax);

	return dif;
}

#ifdef __cplusplus
}
#endif

/*!
 * \brief Read PPM image
 * \param a Double array
 * \param b PPM image file name
 * \return Success: loaded PPM image, Error: a
 */
pRnm read_image (pRnm const &a, const pstring &b) {
	ubyte type, quiet, *dd;
	int k, n, m;
	double *mm;
	PPMimage *image;

	quiet = 1;
	image = loadPPM(b->c_str(), &type, quiet);

	if (!image) {
		std::cerr << DISP_ERROR << " Error loadPPM image " << *b << endl;
		CompileError((string)DISP_ERROR + " Error loadPPM image ");
		return a;
	}

	if (verbosity)
		cout << DISP_INFO << " Size of image : " << image->sizeX << " x " << image->sizeY << " - Type = " << (int)type << endl;

	n = image->sizeX;
	m = image->sizeY;
	a->init(n, m);
	dd = image->data;

	k = 0;
	mm = *a;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			*mm++ = (double)dd[k++] / 256.;

	assert(k == n * m);
	free(image->data);
	free(image);
	return a;
}

/*!
 * \brief Equal operator between array and double array
 * \param a array
 * \param b double array
 * \return a
 */
pRn seta (pRn const &a, const pRnm &b) {
	*a = *b;
	return a;
}

/*!
 * \brief Dynamic load
 */
static void Load_PPM2RNM () {
	cout << " load: ppm2rmn" << endl;

	TheOperators->Add("<-",
	                  new OneOperator2_<KNM<double> *, KNM<double> *, string *>(&read_image)
	                  );
	TheOperators->Add("=",
	                  new OneOperator2_<KN<double> *, KN<double> *, KNM<double> *>(seta)
	                  );
}

LOADFUNC(Load_PPM2RNM)
