/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "Graphics.h"

void GBox()
{
	GStartPath();
	switch(GP->bty) {
	case 'o':
	case 'O':
		GMoveTo(GP->plt[0], GP->plt[2]);
		GLineTo(GP->plt[1], GP->plt[2]);
		GLineTo(GP->plt[1], GP->plt[3]);
		GLineTo(GP->plt[0], GP->plt[3]);
		GLineTo(GP->plt[0], GP->plt[2]);
		break;
	case 'l':
	case 'L':
		GMoveTo(GP->plt[0], GP->plt[3]);
		GLineTo(GP->plt[0], GP->plt[2]);
		GLineTo(GP->plt[1], GP->plt[2]);
		break;
	case '7':
		GMoveTo(GP->plt[0], GP->plt[3]);
		GLineTo(GP->plt[1], GP->plt[3]);
		GLineTo(GP->plt[1], GP->plt[2]);
		break;
	case 'c':
	case 'C':
		GMoveTo(GP->plt[1], GP->plt[2]);
		GLineTo(GP->plt[0], GP->plt[2]);
		GLineTo(GP->plt[0], GP->plt[3]);
		GLineTo(GP->plt[1], GP->plt[3]);
		break;
	default:
		break;
	}
	GEndPath();
}
