=head1 NAME

C<Math::Pari> - Perl interface to PARI.

=head1 SYNOPSIS

  use Math::Pari;
  $a = PARI 2;
  print $a**10000;

or

  use Math::Pari qw(mod);
  $a = mod(3,5);
  print $a**10000;

=head1 DESCRIPTION

This package allows use of most PARI functions as Perl functions. In
what follows we suppose prior knowledge of what PARI is (see
L<ftp://megrez.math.u-bordeaux.fr/pub/pari>).

=head1 EXPORT

=over 4

=item DEFAULT

By default the package exports functions PARI(), PARIcol() and
PARImat() that converts its argument(s) to a 
PARI object. (In fact PARI() is just an alias for C<new Math::Pari>).
The function PARI() accepts following data as its arguments

=over 17

=item One integer

Is converted to a PARI integer.

=item One float

Is converted to a PARI float.

=item One string

Is executed as a PARI expresion (so should not contain whitespace).

=item PARI object

Is passed unchanged.

=item Reference to a Perl array

Each element is converted using the same rules, PARI vector-row with these
elements is returned.

=item Several of above

The same as with a reference to array.

=back

=item Conflicts of rules in PARI()

In deciding what rule of the above to apply the preference is given to
the uppermost choice of available now, if none matches, then the string
rule is used. So C<PARI(1)> returns integer, C<PARI(1.)> returns
float, C<PARI("1")> evaluates "1" as a PARI expression, though all
these data can be converted inside Perl into integer, float or
string. Only what the argument I<is now> is important.

=item PARIcol() and PARImat()

PARIcol() behaves in the same way as PARI() unless given several
arguments. In the latter case it returns a vector-column instead of
vector-row. 

PARImat() constructs a matrix out of the given arguments. It will work
if PARI() will construct a vector of vectors given the same arguments.

=item C<use> with arguments

If arguments are specified in the C<use Math::Pari> directive, the
PARI functions appearing as arguments are exported in the caller
context. In this case the function PARI() and friends is not exported,
so if you need them, you should include them into export list
explicitely, or include C<:DEFAULT> tag. 

The other tags recognized are C<:all>, and number tags C<:4>.  These
tags export functions from the PARI library from the given class
(except for C<:all>, which exports all the classes).  Note that
C<:all> does not include C<:DEFAULT>.

=back

=head1 Available functions

=head2 Directly accessible from Perl

This package supports I<all> the functions from the PARI library with
a signature from a short list. This means that when you update the
PARI library, the newly added function will we available without any change
to this package (provided their signature is in the supported
list). You can reach unsupported functions using
string argument of PARI() function, as in

  3 + PARI('o(x^17)')

A perl script C<parifonc> is provided that lists the functions from
the current release of PARI that are unavailable with the current
release of this glue code. You should specify two arguments to this
script: a path to C<anal.c> file of PARI/GP distribution, and path to
the file  C<Pari.xs> from your truely package.

=head2 Arguments

Arguments to PARI functions are converted to C<long> or PARI type
depending on what type the actual library function requires. No error
checking on arguments is done, so if C<gp> rejects your code since a
particular argument should be of C<type 1> (i.e., a Pari integer),
C<Math::Pari> will silently convert it to C<long>. Each argument is
converted by the rules applicable to PARI().

=head2 Return values

PARI functions return PARI type or a C<long> depending on what the
actual library function returns. 

=head2 Additional functions

Some PARI functions are available in C<gp> (i.e., in C<PARI>
calculator) via infix notation only. In C<Math::Pari> these functions
are available in functional notations too. Some other convenience
functions are also made available.

=over 5

=item Infix, prefix and postfix operations

are available under names

  gneg, gadd, gsub, gmul, gdiv, gdivent, gmod, gpui,
  gle, gge, glt, ggt, geq, gne, gegal, gor, gand,
  gcmp, gcmp0, gcmp1, gcmp_1.

C<gdivent> means euclidean quotient, C<gpui> is power, C<gegal> checks
whether two objects are equal, C<gcmp> is applicable to two real
numbers only, C<gcmp0>, C<gcmp1>, C<gcmp_1> compare with 0, 1 and -1
correspondingly  (see PARI user manual for details).

=item Convenience functions

  pari2iv, pari2nv, pari2num, pari2pv, pari2bool

convert a PARI object to an integer, float, integer/float (whatever is
better), string, and a boolean value correspondingly. Most the time
you do not need these functions due to automatic conversions.

=item Printout functions

  pari_print, pari_pprint, pari_texprint

perform conversions to strings as their PARI counterparts, but do not
print the result.  The difference of pari_print() with pari2pv() is
the number of significant digits they print.)

=item Constant functions

Some mathematical constant appear as function without arguments in
PARI. Corresponding functions are made available under Perl. If you
export them like in

  use Math::Pari wq(:DEFAULT pi i euler);

they can be used as barewords in your program, with the usual
restrictions that sometimes you should disambiguate the call like in

  $negOne = exp(pi() * i);

The parentheses after C<pi> are needed since Perl thinks you want to call
C<pi> with argument C<*i> otherwise (since C<*i> is a legal Perl
expression, and function calls are greedy on the right).

=back

=head1 PARI objects

Functions from PARI library take as arguments and/or return objects of
type C<GEN> (in C<C> notations). In Perl these data are encapsulated
into special kind of Perl variables: PARI objects. You can check for a
variable C<$obj> to be a PARI object using

  ref $obj eq 'Math::Pari';

Most the time you do not need this due to automatic conversions.

=head1 PARI polynomials and Perl barewords

Some bareletters denote Perl operators, like C<q>, C<x>, C<y>,
C<s>. This can lead to errors in Perl parsing your expression. Say, while

  print sin(tan(x))-tan(sin(x))-asin(atan(x))+atan(asin(x));

parses OK (after C<use Math::Pari qw(sin tan asin atan)>),

  print sin(tan(y))-tan(sin(y))-asin(atan(y))+atan(asin(y));

does not. You should avoid lower-case barewords used as PARI variables.

=head1 Overloading and automatic conversion

Whenever an arithmetic operation includes a PARI object the other
arguments are converted to a PARI type and the corresponding PARI
library functions is used to implement the operation. Numeric
comparison operations use C<gcmp> and friends, string comparisons compare in
lexicographical order using C<lex>. Currently the following arithmetic
operations are overloaded:

  unary -, +, -, *, /, %, **, abs, cos, sin, exp, log, sqrt.

Whenever a PARI object appears in a situation that requires integer,
numeric, boolean or string data, it is converted to the corresponding
type. Boolean conversion is subject to usual PARI pitfalls related to
imprecise zeros (see documentation of C<gcmp0> in PARI reference).

=head1 PREREQUISITES

=head2 Perl

To compile the extension you need to have at least C<MakeMaker>
version 3.7. In the versions of perl earlier than 5.001 negative
constants were converted to floats, so to use PARI operations that do
different things on integers and floats you would like to use a more
recent version of perl. Overloading was buggy in 5.000.

=head2 PARI

You need at least version 1.39 of PARI. (See
L<ftp://megrez.math.u-bordeaux.fr/pub/pari>.)

=head1 Perl vs. PARI: different syntax

Note that the PARI notations should be used in string arguments to
PARI() function, while Perl notations should be used otherwise.

=over 4

=item C<^>

Power is denoted by C<**> in Perl.

=item C<\> and C<\/>

There are no such operators in Perl, use the word forms
C<gdivent(x,y)> and C<gdivround(x,y)> instead.

=item C<~>

There is no postfix C<~> Perl operator.  Use trans() instead.

=item C<!>

There is no postfix C<!> Perl operator.  Use fact() instead (note that
the semantics of C<fact(x)> and C<x!> is different in PARI, one is
real, another integer).

Currently there is no direct way to get integer factorial in Math::Pari.

=item doubles

Doubles in Perl are of precision approximately 15 digits.  When you
use them as arguments to PARI functions, they are converted to PARI
real variables, and due to intermediate 15-decimal-to-binary
conversion of Perl variables the result may be different than with the
PARI many-decimal-to-binary conversion.  Say, C<PARI(0.01)> and
C<PARI('0.01')> differ at 19-th place, as

  setprecision(38); print pari_print(0.01), "\n", pari_print('0.01'), "\n";

shows.

=back

=head1 ENVIRONMENT

No environment variables are used.

=head1 BUGS

=over 5

=item *

Not all the PARI functions are directly available.

=item *

Many others...

=back

=head1 AUTHOR

Ilya Zakharevich, I<ilya@math.mps.ohio-state.edu>

=cut

# $Id: Pari.pm,v 1.3 1994/11/25 23:40:52 ilya Exp ilya $
package Math::Pari;

require Exporter;
require DynaLoader;

@ISA = qw(Exporter AutoLoader DynaLoader);

# Items to export into callers namespace by default
# (move infrequently used names to @EXPORT_OK below)

@EXPORT = qw(
PARI PARIcol PARImat
);

# Other items we are prepared to export if requested (may be extended during
# ->import. )
@EXPORT_OK = qw(
  sv2pari sv2parimat pari2iv pari2nv pari2num pari2pv pari2bool loadPari _bool
  listPari pari_print pari_pprint pari_texprint
);

use subs qw(
   _gneg
   _gadd
   _gsub
   _gmul
   _gdiv
   _gmod
   _gpui
   _gle
   _gge
   _glt
   _ggt
   _geq
   _gne
   _gcmp
   _lex
   _2bool
   pari2pv
   pari2num
   _abs
   _cos
   _sin
   _exp
   _log
   _sqrt
);				# For overloading

use overload qw(
   neg _gneg
   + _gadd
   - _gsub
   * _gmul
   / _gdiv
   % _gmod
   ** _gpui
   <= _gle
   >= _gge
   < _glt
   > _ggt
   == _geq
   != _gne
   <=> _gcmp
   cmp _lex
   bool _2bool
   "" pari2pv
   +0 pari2num
   abs _abs
   cos _cos
   sin _sin
   exp _exp
   log _log
   sqrt _sqrt
);

## Was needed with the old way of overloading
## foreach (values %OVERLOAD) {$interchange{$_}=$' if /^_/}
# Against a bug in autoloading
foreach (keys %OVERLOAD) {$OVERLOAD{$_}= \&{ $OVERLOAD{$_} }}

sub AUTOLOAD {
  $AUTOLOAD =~ /^(?:Math::Pari::)?/;
# #   my $name1=$';
# #   my $name=$name1;
# #   if ($interchange{$name1}) {
# #     $name = $interchange{$name1};
# #     eval <<EOE;
# #       sub $name1 {
# # 	\$_[2]? $name(\$_[1],\$_[0]):  $name(\$_[0],\$_[1]);
# #       }
# # EOE
# #   }
  my $cv=loadPari($');
  
#  goto &$cv;
#  goto &$AUTOLOAD;
#  &$cv;
  &{$'}(@_);
#  &$AUTOLOAD;
}

# Needed this guy to circumwent autoloading while no XS definition

#### sub DESTROY {}


# sub AUTOLOAD {
#     if ((caller(0))[4]) {
# 	$AutoLoader::AUTOLOAD = $AUTOLOAD;
# 	goto &AutoLoader::AUTOLOAD;
#     }
#     local($constname);
#     ($constname = $AUTOLOAD) =~ s/.*:://;
#     $val = constant($constname, @_ ? $_[0] : 0);
#     if ($! != 0) {
# 	if ($! =~ /Invalid/) {
# 	    $AutoLoader::AUTOLOAD = $AUTOLOAD;
# 	    goto &AutoLoader::AUTOLOAD;
# 	}
# 	else {
# 	    ($pack,$file,$line) = caller;
# 	    die "Your vendor has not defined Math::Pari macro $constname, used at $file line $line.
# ";
# 	}
#     }
#     eval "sub $AUTOLOAD { $val }";
#     goto &$AUTOLOAD;
# }

$initmem = $initmem || 4000000;		# How much memory for the stack
$initprimes = $initprimes || 500000;	# Calculate primes up to this number

bootstrap Math::Pari;

# Preloaded methods go here.  Autoload methods go after __END__, and are
# processed by the autosplit program.

sub new {
  shift;
  if (@_>1) {my(@t)=@_;return sv2pari(\@t)}
  return sv2pari(shift);
}

###sub PARI {new Math::Pari @_}

sub import {
  my $p=shift;
  @_ = map {
    if (/^:(?!DEFAULT)(.*)/) {
      my $tag = $1;
      $tag = -1 if ($tag eq 'all');
      listPari($tag);
    } else {
      ($_);
    }
  } @_;
  
  push( @EXPORT_OK,grep( (eval {loadPari($_)}, !$@) , @_) );
  #warn "EXPORT_OK: `@EXPORT_OK', import called with: `@_'\n";
  #Exporter::import $p @_;
  &Exporter::export($p,(caller(0))[0],@_);
}

1;
__END__
