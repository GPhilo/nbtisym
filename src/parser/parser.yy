#line 11148 "./doc/bison.texi"
%skeleton "lalr1.cc" /* -*- C++ -*- */
%require "3.0"
%defines
%define parser_class_name {calcxx_parser}
#line 11164 "./doc/bison.texi"
%define api.token.constructor
%define api.value.type variant
%define parse.assert
#line 11181 "./doc/bison.texi"
%code requires
{
# include <string>
class calcxx_driver;
}
#line 11195 "./doc/bison.texi"
// The parsing context.
%param { calcxx_driver& driver }
#line 11207 "./doc/bison.texi"
%locations
%initial-action
{
  // Initialize the initial location.
  @$.begin.filename = @$.end.filename = &driver.file;
};
#line 11222 "./doc/bison.texi"
%define parse.trace
%define parse.error verbose
#line 11233 "./doc/bison.texi"
%code
{
# include "driver.hh"
}
#line 11249 "./doc/bison.texi"
%define api.token.prefix {TOK_}
%token
  END  0  "end of file"
  COLON
  ID_FIELD
  DC_FIELD
  VDD_FIELD
  FREQ_FIELD
;
#line 11269 "./doc/bison.texi"
%token <std::string> IDENTIFIER "identifier"
%token <float> NUMBER "number"
%type <std::string> ID
%type <float> VDD
%type <float> FREQ
%type <float> DC
#line 11282 "./doc/bison.texi"
%printer { yyoutput << $$; } <*>;
#line 11291 "./doc/bison.texi"
%%
%start states;

states	: state				{}
		| states state		{}
		;

state	: ID VDD FREQ DC	{ driver.states[$1] = ProcessorState($1, $2, $3, $4); }
		;

ID		: ID_FIELD COLON IDENTIFIER		{ $$ = $3; }
		;

VDD		: VDD_FIELD COLON NUMBER		{ $$ = $3; }
		;

FREQ		: FREQ_FIELD COLON NUMBER		{ $$ = $3; }
		;

DC		: DC_FIELD COLON NUMBER		{ $$ = $3; }
		;
%%
#line 11321 "./doc/bison.texi"
void
yy::calcxx_parser::error (const location_type& l,
                          const std::string& m)
{
  driver.error (l, m);
}
