/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     STRING = 258,
     ID = 259,
     NUM = 260,
     LBRACK = 261,
     RBRACK = 262,
     ACCELERATOR = 263,
     AREALIGHTSOURCE = 264,
     ATTRIBUTEBEGIN = 265,
     ATTRIBUTEEND = 266,
     CAMERA = 267,
     CONCATTRANSFORM = 268,
     COORDINATESYSTEM = 269,
     COORDSYSTRANSFORM = 270,
     FILM = 271,
     IDENTITY = 272,
     LIGHTSOURCE = 273,
     LOOKAT = 274,
     MATERIAL = 275,
     OBJECTBEGIN = 276,
     OBJECTEND = 277,
     OBJECTINSTANCE = 278,
     PIXELFILTER = 279,
     REVERSEORIENTATION = 280,
     ROTATE = 281,
     SAMPLER = 282,
     SCALE = 283,
     SEARCHPATH = 284,
     SHAPE = 285,
     SURFACEINTEGRATOR = 286,
     TEXTURE = 287,
     TRANSFORMBEGIN = 288,
     TRANSFORMEND = 289,
     TRANSFORM = 290,
     TRANSLATE = 291,
     VOLUME = 292,
     VOLUMEINTEGRATOR = 293,
     WORLDBEGIN = 294,
     WORLDEND = 295,
     HIGH_PRECEDENCE = 296
   };
#endif
/* Tokens.  */
#define STRING 258
#define ID 259
#define NUM 260
#define LBRACK 261
#define RBRACK 262
#define ACCELERATOR 263
#define AREALIGHTSOURCE 264
#define ATTRIBUTEBEGIN 265
#define ATTRIBUTEEND 266
#define CAMERA 267
#define CONCATTRANSFORM 268
#define COORDINATESYSTEM 269
#define COORDSYSTRANSFORM 270
#define FILM 271
#define IDENTITY 272
#define LIGHTSOURCE 273
#define LOOKAT 274
#define MATERIAL 275
#define OBJECTBEGIN 276
#define OBJECTEND 277
#define OBJECTINSTANCE 278
#define PIXELFILTER 279
#define REVERSEORIENTATION 280
#define ROTATE 281
#define SAMPLER 282
#define SCALE 283
#define SEARCHPATH 284
#define SHAPE 285
#define SURFACEINTEGRATOR 286
#define TEXTURE 287
#define TRANSFORMBEGIN 288
#define TRANSFORMEND 289
#define TRANSFORM 290
#define TRANSLATE 291
#define VOLUME 292
#define VOLUMEINTEGRATOR 293
#define WORLDBEGIN 294
#define WORLDEND 295
#define HIGH_PRECEDENCE 296




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 135 "core/pbrtparse.y"
{
char string[1024];
float num;
ParamArray *ribarray;
}
/* Line 1529 of yacc.c.  */
#line 137 "core/pbrtparse.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

