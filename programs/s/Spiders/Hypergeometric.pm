#!/usr/bin/perl

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: Spiders::Hypergeometric

package Spiders::Hypergeometric;
        
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new {
	my ($class,%arg)=@_;
    my $self = {
        _n => undef,
		_k => undef,
		_r =>undef,
    };
    bless $self, $class;
	return $self;
}

sub Hypergeometric {
	my ($self, $n, $k, $r, $x) = @_;
	return $k == 0 if $r == 0;
	return exp(choose($r, $x) + choose($n - $r, $k - $x) - choose($n, $k));
}

sub choose {
	my ($n, $k) = @_;
	my ($result, $j) = (0, 1);
	return 0 if $k > $n || $k < 0;
	$k = ($n - $k) if ($n - $k) < $k;
	while ($j <= $k) {
		$result += log($n--);
		$result -= log($j++);
	}
	return $result;
}

=head
sub logfact {
   return gammln(shift(@_) + 1.0);
}

sub hypergeom {

   my ($self,$m, $n, $N, $i) = @_;

   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}

sub gammln {
  my $xx = shift;
  my @cof = (76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.12086509738661e-2, -0.5395239384953e-5);
  my $y = my $x = $xx;
  my $tmp = $x + 5.5;
  $tmp -= ($x + .5) * log($tmp);
  my $ser = 1.000000000190015;
  for my $j (0..5) {
     $ser += $cof[$j]/++$y;
  }
  -$tmp + log(2.5066282746310005*$ser/$x);
}
=cut

1;
