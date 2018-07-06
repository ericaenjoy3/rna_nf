#!/usr/bin/env perl
use strict;
#Perl script to sort intervals in STDIN based on given chromosome order
#By Minghui Wang

my %chr_order;
my $n=0;
open(FI,"$ARGV[0]") || die "can't open $ARGV[0]\n";

while(<FI>){
	if(/^\@SQ\tSN:(\S+)\tLN:\d+/){
		$chr_order{$1}=$n++;
	}
}
close(FI);

my %interval;
while(<STDIN>){
	my @arr=split(/\t/,$_);
	$interval{$arr[0]}{$arr[1]}{$arr[2]}=$_;
}

for my $chr(sort {$chr_order{$a} <=> $chr_order{$b}} keys %interval){
	for my $s1(sort {$a <=> $b} keys %{$interval{$chr}}){
		for my $s2(sort {$a <=> $b} keys %{$interval{$chr}{$s1}}){
			print STDOUT $interval{$chr}{$s1}{$s2};
		}
	}
}
