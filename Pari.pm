=head1 NAME

C<Math::Pari> - Perl interface to PARI.

=head1 SYNOPSIS

  use Math::Pari;
  $a = PARI 2;
  print $a**10000;

or

  use Math::Pari qw(Mod);
  $a = Mod(3,5);
  print $a**10000;

=head1 DESCRIPTION

This package is a Perl interface to famous library PARI for
numerical/scientific/number-theoretic calculations.  It allows use of
most PARI functions as Perl functions, and (almost) seamless merging
of PARI and Perl data. In what follows we suppose prior knowledge of
what PARI is (see L<ftp://megrez.math.u-bordeaux.fr/pub/pari>, or
L<Math::libPARI>).

=head1 EXPORTed functions

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
the uppermost choice of those available now.  If none matches, then the string
rule is used. So C<PARI(1)> returns integer, C<PARI(1.)> returns
float, C<PARI("1")> evaluates "1" as a PARI expression.

Note that for Perl these data are synonimous, since Perl freely
converts between integers, float and strings.  However, to PARI() only
what the argument I<is now> is important.

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
explicitely, or include C<:DEFAULT> tag:

  use Math::Pari qw(factorint PARI);
  use Math::Pari qw(:DEFAULT factorint);

or simply do it in two steps

  use Math::Pari;  
  use Math::Pari 'factorint';

The other tags recognized are C<:PARI>, C<:all>, C<prec=NUMBER>,
number tags, like C<:4>, and section names tags.  The number tags
export functions from the PARI library from the given class (except
for C<:PARI>, which exports all the classes).  Tag C<:all> exports all
the exportable symbols and C<:PARI>.

Giving C<?> command to C<gp> B<PARI> calculator lists the following classes: 

  1: Standard monadic or dyadic OPERATORS
  2: CONVERSIONS and similar elementary functions
  3: TRANSCENDENTAL functions
  4: NUMBER THEORETICAL functions
  5: Functions related to ELLIPTIC CURVES
  6: Functions related to general NUMBER FIELDS
  7: POLYNOMIALS and power series
  8: Vectors, matrices, LINEAR ALGEBRA and sets
  9: SUMS, products, integrals and similar functions
  10: GRAPHIC functions
  11: PROGRAMMING under GP

One can use section names instead of number tags.  Recognized names are

  :standard :conversions :transcendental :number :elliptic
  :fields :polynomials :vectors :sums :graphic :programming

One can get a list of all C<Math::Pari> accessible functions, or
functions from the given section using listPari() function.

Starting from version 5.005 of Perl, three additional tags are
supported: C<:int>, C<:float>, C<:hex>.  If used, all the
integer/float/hex-or-octal literals in Perl will be automatically
converted to became PARI objects.  Say

  use Math::Pari ':int';
  print 2**1000;

is equivalent to

  print PARI(2)**PARI(1000);

=back

=head1 Available functions

=head2 Directly accessible from Perl

This package supports I<all> the functions from the PARI library with
a I<signature> which can be recognized by Math::Pari.  This means that
when you update the PARI library, the newly added function will we
available without any change to this package (provided their signature
is supported).  You can reach unsupported functions using string
argument of PARI() function, as in

  3 + PARI('O(x^17)')

(or some special wrapper functions, like C<O(variable,power)>).
A perl script C<parifonc> is provided that lists the functions from
the current release of PARI that are unavailable with the current
release of this glue code (not checked with PARI 2.0).

The following functions are specific to GP calculator, thus are not
available to Math::Pari in any way:

  allocatemem default error extern input print print1 printp printp1
  printtex quit read system whatnow write write1 writetex

whatnow() function is useless, since Math::Pari does not support the
"compatibility" mode (with older PARI library).  The functionality of
print(), write() and variants is available via automatic string
translation, and pari_print() function and its variants.

allocatemem() and default() are the only two important functions with
functionality not supported by the current interface.  Note however,
that two most important default() actions are supported by
setprecision() and setseriesprecision() functions.

=head2 Arguments

Arguments to PARI functions are converted to C<long> or PARI type
depending on what type the actual library function requires. No error
checking on arguments is done, so if C<gp> rejects your code since a
particular argument should be of C<type T_INT> (i.e., a Pari integer),
the corresponding function of C<Math::Pari> will silently convert the
argument to C<long>. Each argument is converted by the rules
applicable to PARI.

=head2 Return values

PARI functions return PARI type or a Perl's integer depending on what
the actual library function returns.

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
correspondingly (see PARI user manual for details, or
L<Math::libPARI>).  Note that all these functions are more readily
available via operator overloading, so instead of

  gadd($x, gneg($y))

one can write

  $x+(-$y)

(as far as overloading may be triggered, so we assume that $x or $y is
of PARI type already).

=item Conversion functions

  pari2iv, pari2nv, pari2num, pari2pv, pari2bool

convert a PARI object to an integer, float, integer/float (whatever is
better), string, and a boolean value correspondingly. Most the time
you do not need these functions due to automatic conversions.

=item Printout functions

  pari_print, pari_pprint, pari_texprint

perform conversions to strings as their PARI counterparts, but do not
print the result.  The difference of pari_print() with pari2pv() is
the number of significant digits they print, and whitespace in the
output.)

=item Constant functions

Some mathematical constant appear as function without arguments in
PARI.  Perl has a facility to have similar functions. If you
export them like in

  use Math::Pari qw(:DEFAULT Pi I Euler);

they can be used as barewords in your program:

  $x = Pi ** Euler;

=item Low-level functions

For convenience of low-level PARI programmers some low-level functions
are made available as well (they are not exported):

  typ(x) changevalue(name,newvalue)

=item Uncompatible functions

  O

Since implementing C<O(7**6)> would be very tedious, we provide a
two-argument form C<O(7,6)> instead.  Note that with polynomials there
is no problem like this one, both C<O($x,6)> and C<O($x**6)> work.

  ifact(n)

integer factorial functions, available from C<gp> as C<n!>.

=back

=head2 Looping functions

PARI has a big collection of functions which loops over some set.
Such a function takes two I<special> arguments: loop variable, and the
code to execute in the loop.

The code can be either a string (which contains PARI code to execute -
thus should not contain whitespace), or a Perl code reference.  The
loop variable can be a string giving the name of PARI variable (as in

  fordiv(28, 'j', 'a=a+j+j^2');

or

  $j= 'j';
  fordiv(28, $j, 'a=a+j+j^2');

), or a Perl variable containing a PARI variable (as in

  $j = PARI 'j';
  fordiv(28, $j, sub { $a += $j + $j**2 });

).  

If the loop variable is not of these two types, then an appropriate
name will be autogenerated.  Note that since you have no control over
this name, you will not be able to use this variable from your PARI
code, say

  $j = 7.8;
  fordiv(28, $j, 'a=a+j+j^2');

will not (obviously) expand C<j> to mirror $j (unless you set up C<j>
to mirror $j explicitely, see L<"Accessing Perl functions from PARI code">).

B<Useless musing alert! Do not read the rest of this section!>

In fact a very hairy type of access is also supported.  Note that the
following code will not do what you expect

  $x = 0;
  $j = PARI 'j';
  fordiv(28, 'j', sub { $x += $j } );

since C<fordiv> will I<localize> C<j> inside the loop, so $j will
still reference the old value, which is an independent variable, not
the index of the loop.  The simplest workaround is not to use the
above syntax (i.e., not mixing literal loop variable with Perl loop
code, just using $j as the second argument to C<fordiv> is enough):

  $x = 0;
  $j = PARI 'j';
  fordiv(28, $j, sub { $x += $j } );

However, if absolutely required, one can make a I<delayed> variable $j
which will always reference the same thing C<j> references in PARI
I<now> by using C<PARIvar> constructor

  $x = 0;
  $j = PARIvar 'j';
  fordiv(28, 'j', sub { $x += $j } );

This problem is related to

  $ref = \$_;			# $$ref is going to be old value even after
				# localizing $_ in Perl's grep/map

not accessing localized values of $_ in the plain Perl.


=head2 Accessing Perl functions from PARI code

This is possible.  Just use the same name for the function:

  sub counter { $i += shift; }
  $i = 145;
  PARI 'k=5' ;
  fordiv(28, 'j', 'k=k+counter(j)');
  print PARI('k'), "\n";

prints 

   984

Note that if the subroutine takes a variable number of arguments, C<@>
in the prototype (or a missing prototype) counts as 6 optional
arguments are supported.  If called from PARI with fewer arguments
optional arguments will be set to integer PARI 0.

Note also that no direct import of Perl variables is available yet
(but you can write a function wrapper for this):

  sub getv () {$v}

There is an unsupported function for explicitely importing Perl
functions into Pari, possibly with a different name, and possibly with
explicitely specifying number of arguments.

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

may parse OK (after C<use Math::Pari qw(sin tan asin atan)>),

  print sin(tan(y))-tan(sin(y))-asin(atan(y))+atan(asin(y));

does not. You should avoid lower-case barewords used as PARI
variables, say, do

  $y = PARI('y');
  print sin(tan($y))-tan(sin($y))-asin(atan($y))+atan(asin($y));

to get

  -1/18*y^9+26/4725*y^11-41/1296*y^13+328721/16372125*y^15+O(y^16)

Well, the best advice: do not use barewords anywhere in your program!

=head1 Overloading and automatic conversion

Whenever an arithmetic operation includes a PARI object the other
arguments are converted to a PARI type and the corresponding PARI
library functions is used to implement the operation. Numeric
comparison operations use C<gcmp> and friends, string comparisons compare in
lexicographical order using C<lex>. Currently the following arithmetic
operations are overloaded:

  unary -
  + - * / % ** abs cos sin exp log sqrt
  << >>
  <= == => <  >  != <=> 
  le eq ge lt gt ne cmp

Whenever a PARI object appears in a situation that requires integer,
numeric, boolean or string data, it is converted to the corresponding
type. Boolean conversion is subject to usual PARI pitfalls related to
imprecise zeros (see documentation of C<gcmp0> in PARI reference).

Note that a check for equality is subject to same pitfalls as in PARI
due to imprecise values.  PARI may also refuse to compare data of
different types for equality if it thinks this may lead to
counterintuitive results.

Note also that numeric ordering is not defined for some types of PARI
objects.  For string comparison operations we use PARI-lexicographical
ordering.

=head1 PREREQUISITES

=head2 Perl

In the versions of perl earlier than 5.003 overloading used a
different interface, so you may need to convert C<use overload> line
to C<%OVERLOAD>, or, better, upgrade.

=head2 PARI

Starting from version 2.0, this module comes without a PARI library included.

For the source of PARI library see
L<ftp://megrez.math.u-bordeaux.fr/pub/pari>.

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

There is no postfix C<~> Perl operator.  Use mattranspose() instead.

=item C<'>

There is no postfix C<'> Perl operator.  Use deriv() instead.

=item C<!>

There is no postfix C<!> Perl operator.  Use factorial()/ifact() instead
(returning a real or an integer correspondingly).

=item big integers

Currently Perl will convert big I<literal> integers to doubles if they
could not be put into B<C> 32-bit signed integers.  If you want to
input such an integer, use PARI('12345678901234567890').

Starting from version 5.005 of Perl, if the tag C<:int> is used, all
the integer literals in Perl will be automatically converted to became
PARI objects.  Say

  use Math::Pari ':int';
  print 2**1000;

is equivalent to

  print PARI(2)**PARI(1000);

=item doubles

Doubles in Perl are of precision approximately 15 digits.  When you
use them as arguments to PARI functions, they are converted to PARI
real variables, and due to intermediate 15-decimal-to-binary
conversion of Perl variables the result may be different than with the
PARI many-decimal-to-binary conversion.  Say, C<PARI(0.01)> and
C<PARI('0.01')> differ at 19-th place, as

  setprecision(38);
  print pari_print(0.01),   "\n",
        pari_print('0.01'), "\n";

shows.

Note that setprecision() changes the output format of pari_print() and
friends, as well as internal precision.  The generic PARI===>string
conversion does not take into account the output format, thus

  setprecision(38);
  print PARI(0.01),       "\n",
        PARI('0.01'),     "\n",
        pari_print(0.01), "\n";

will print all the lines with different number of digits after the
point: the first one with 22, since the double 0.01 was converted to a
low-precision PARI object, the second one with 41, since internal form
for precision 38 requires that many digits for representation, and the
last one with 39 to have 38 significant digits.

Starting from version 5.005 of Perl, if the tag C<:float> is used, all
the float literals in Perl will be automatically converted to became
PARI objects.  Say

  use Math::Pari ':float';
  print atan(1.);

is equivalent to

  print atan(PARI('1.'));

=item array base

Arrays are 1-based in PARI, are 0-based in Perl.  So while array
access is possible in Perl, you need to use different indices:

  $nf = PARI 'nf';	# number field
  $a = PARI('nf[7]');
  $b = $nf->[6];

Now $a nd $b contain the same value.

=item matrices

Note that C<PARImat([[...],...,[...])> constructor creates a matrix
with specified columns, while PARI's C<[1,2,3;4,5,6]> constructor
creates a matrix with specified rows.  Use a convenience function
PARImat_tr() which will transpose a matrix created by PARImat() to use
the same order of elements as in PARI.

=item builtin perl functions

Some PARI functions, like C<length> and C<eval>, are Perl
(semi-)reserved words.  To reach these functions, one should either
import them: 

  use Math::Pari qw(length eval);

or call them with prefix (like C<&length>) or the full name (like
C<Math::Pari::length>).

=item string($w, $text)

will not coerce a number into C<8+3> format.  You need to use
C<sprintf "%8.3f", $text> explicitely.

=back

=head1 High-resolution graphics

If you have Term::Gnuplot installed, you may use high-resolution
graphic primitives of B<PARI>.  Before the usage you need to establish
a link between Math::Pari and Term::Gnuplot by calling link_gnuplot().
You can change the output filehandle by calling set_plot_fh(), and
output terminal by calling plotterm(), as in

    use Math::Pari qw(:graphic asin);

    open FH, '>out.tex' or die;
    link_gnuplot();
    set_plot_fh(\*FH);
    plotterm('emtex');
    ploth($x, .5, .999, sub {asin $x});
    close FH or die;

=head1 libPARI documentation

libPARI documentation is included, see L<Math::libPARI>.  It is converted
from Chapter 3 of B<PARI/GP> documentation by F<paridoc_to_pod> script.

=head1 ENVIRONMENT

No environment variables are used.

=head1 BUGS

=over 5

=item *

A few of PARI functions are available indirectly only.

=item F<t/*_will_fail.t>

These tests expose several bugs.

=back

=head1 AUTHOR

Ilya Zakharevich, I<ilya@math.ohio-state.edu>

=cut

# $Id: Pari.pm,v 1.3 1994/11/25 23:40:52 ilya Exp ilya $
package Math::Pari::Arr;

#sub TIEARRAY { $_[0] }
sub STORE { die "Storing into array elements unsupported" }

package Math::Pari;

require Exporter;
require DynaLoader;
#use autouse Carp => 'croak';

@ISA = qw(Exporter DynaLoader);
@Math::Pari::Ep::ISA = qw(Math::Pari);

# Items to export into callers namespace by default
# (move infrequently used names to @EXPORT_OK below)

@EXPORT = qw(
PARI PARIcol PARImat PARIvar PARImat_tr
);

# Other items we are prepared to export if requested (may be extended during
# ->import. )
@EXPORT_OK = qw(
  sv2pari sv2parimat pari2iv pari2nv pari2num pari2pv pari2bool loadPari _bool
  listPari pari_print pari_pprint pari_texprint O ifact gdivent gdivround
  changevalue set_plot_fh link_gnuplot setprecision setseriesprecision
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

sub _shiftl {
  my ($left,$right) = (shift,shift);
  ($left,$right) = ($right, $left) if shift;
  $left * 2**$right;
}

sub _shiftr {
  my ($left,$right) = (shift,shift);
  ($left,$right) = ($right, $left) if shift;
  floor($left / 2**$right);
}

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
   <<  _shiftl
   >>  _shiftr
);

sub AUTOLOAD {
  $AUTOLOAD =~ /^(?:Math::Pari::)?(.*)/;
  my $cv = loadPari($1);
  
#  goto &$cv;
#  goto &$AUTOLOAD;
#  &$cv;
  &$1;
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

$VERSION = '2.001403';

bootstrap Math::Pari;

# Preloaded methods go here.  Autoload methods go after __END__, and are
# processed by the autosplit program.

sub new {
  shift;
  if (@_>1) {my(@t)=@_;return sv2pari(\@t)}
  return sv2pari(shift);
}

###sub PARI {new Math::Pari @_}

%names    = qw(
	       1 standard
	       2 conversions
	       3 transcendental
	       4 number
	       5 elliptic
	       6 fields
	       7 polynomials
	       8 vectors
	       9 sums
	       10 graphic
	       11 programming
	      );
@sections{values %names} = keys %names;

@converted{split /,\s+/, qq(buchimag, buchreal,
    buchgen, buchgenforcefu, buchgenfu, buchinit, buchinitforcefu, buchinitfu,
    string, addhelp, kill)} = (1) x 100;

# Highly unfinished
sub _cvt { PARI(shift) }
sub _hex_cvt {
  my $in = shift;
  my $mult = PARI(1);
  my $ret = 0;
  my $shift = 1<<(4*7);

  $in =~ s/^0(x)?// or die;
  my $hex = $1;
  $shift = 1<<(3*7) unless $hex;
  while ($in =~ s/([a-fA-F\d]{1,7})$//) {
    my $part = $hex ? hex $1 : oct $1;
    
    $ret += $part * $mult;
    $mult *= $shift;
  }
  die "Cannot hex '$in'" if length $in;
  return $ret;
}
%overloaded_const = ( 'int' => \&_cvt, float => \&_cvt, 'hex' => \&_hex_cvt);
%overloaded_const_word
  = ( 'int' => 'integer', float => 'float', 'hex' => 'binary');

sub import {
  my $p=shift;
  my @consts;			# Need to do it outside any block!
  @_ = map {
    if (/^:(?!DEFAULT)(.*)/) {
      my $tag = $1;
      my $sect = $tag;
      my @pre;
      $tag = -1, @pre = (@EXPORT_OK,@EXPORT) if ($tag eq 'all');
      $tag = -1 if ($tag eq 'PARI');
      $tag = $sections{$tag} if $tag !~ /^-?\d+$/ and exists $sections{$tag};
      push @pre, 'link_gnuplot', 'set_plot_fh' if $tag eq 10;
      if ($tag =~ /^prec=(\d+)$/) {
	setprecision($1);
	();
      } elsif ($tag =~ /^(int|hex|float)$/) {
	die "Overloaded constants are not supported in this version of Perl"
	  if $] < 5.004_69;
	push @const, $overloaded_const_word{$tag} => $overloaded_const{$tag};
	# print "Constants: $overloaded_const_word{$tag} => $overloaded_const{$tag} \n";
	();
      } elsif (defined $tag and $tag =~ /^-?\d+$/) {
	(@pre, listPari($tag));
      } else {
	die "Unknown section '$sect' specified";
      }
    } else {
      ($_);
    }
  } @_;

  overload::constant(splice @const, 0, 2) while @const;  

  # print "EXPORT_OK: @EXPORT_OK\n";
  push @EXPORT_OK,
      grep( ($_ ne ':DEFAULT' 
	     and not $export_ok{$_}
	     and (eval {loadPari($_), 1} or warn $@), !$@) ,
	    @_);
  # print "EXPORT_OK: @EXPORT_OK\n";
  &Exporter::export($p,(caller(0))[0],@_);
}

sub O ($;$) {
  return PARI("O($_[0]^$_[1])") if @_ == 2;
  return PARI("O($_[0])") if typ($_[0]) == 10; # Poly
  Carp::croak("O(number**power) not implemented, use O(number,power) instead");
}

sub PARImat_tr {mattranspose(PARImat(@_))}
sub string ($$) {
  PARI (qq'string($_[0],"$_[1]")');
}

sub installPerlFunction {my @a=@_; $a[0] = \&{$a[0]}; installPerlFunctionCV(@a)}

my $name;

for $name (keys %converted) {
  push @EXPORT_OK, $name;
  next if defined &$name;
  # string needs to format numbers to 8.3...
  if ($name eq 'addhelp' or $name eq 'string') {
    *$name = sub {PARI("$name($_[0],'$_[1]')")}
  } else {
    *$name = sub {local $"=',';PARI("$name(@_)")}
  }
}

@export_ok{@EXPORT_OK,@EXPORT} = (1) x (@EXPORT_OK + @EXPORT);

sub link_gnuplot {
  eval 'use Term::Gnuplot 0.4; 1' or die;
  set_gnuterm(Term::Gnuplot::change_term_address(), 
	      Term::Gnuplot::term_tbl_address());
}

sub set_plot_fh {
  eval 'use Term::Gnuplot 0.4; 1' or die;
  Term::Gnuplot::set_gnuplot_fh(@_);
}

1;
__END__
