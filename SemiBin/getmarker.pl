#script from Maxbin2 https://sourceforge.net/projects/maxbin2/
#!/usr/bin/perl -w
use strict;

my $COV_CUTOFF = 0.4;

gethmmmarker($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

# void gethmmmarker(input hmm output, fasta file, min seq length, output, is_harder);
sub gethmmmarker
{

	my $hmm = $_[0];
	my $fasta_f = $_[1];
	my $min_seq_len = $_[2];
	my $out = $_[3];
	my $localine;
	my @arr;
	my $i;
	my $j;
	my $inlinehmm;
	my $current = "";
	my %lenhash;
	my %tmphash;
	my %queryhash;
	my $querycount = 0;
	my @queryseq;
	my @queryseqnum;
	my @querylen;
	my $tmp1;

	# Read fasta file
	open(LOCALFILE, "<$fasta_f");
	while(defined($localine = <LOCALFILE>))
	{
		chomp($localine);
		if ($localine =~ /^>/)
		{
			$current = trim(substr($localine, 1));
			$lenhash{$current} = 0;
		}
		elsif ($localine ne "")
		{
			$lenhash{$current} = $lenhash{$current} + length($localine);
		}
	}
	close(LOCALFILE);

	# Read HMM output
	$current = "";
	open(LOCALFILE, "<$hmm") || die "Cannot open hmm output file $hmm\n";
	while(defined($localine = <LOCALFILE>))
	{	
		chomp($localine);
		if (index($localine, "#") == 0)
		{
			next;
		}
		@arr = split(/[ ]+/, $localine);
		#$inlinehmm = checkMarker($arr[4]);
		$inlinehmm = checkMarker($arr[3]);
		#print "$arr[0], $inlinehmm, $arr[5], $arr[15], $arr[16]\n";
		#if ($arr[0] =~ /([A-Za-z0-9._:;'"`=~\:\!\@\#\$\%\^\&\*\(\)\{\}\[\]\\\/\?\<\>\-\|]+)_[0-9]+$/)
		if ($arr[0] =~ /([A-Za-z0-9._:;'"`=~\:\!\@\#\$\%\^\&\*\(\)\{\}\[\]\\\/\?\<\>\-\|]+)_[0-9]+_[0-9]+_[+\-]$/)
		{
			if (exists $lenhash{$1} && $lenhash{$1} >= $min_seq_len)
			{
				if ($current ne $inlinehmm)
				{
					if ($current ne "")
					{
						if ((scalar(keys %tmphash)) > 0)
						{
							# Flush tmphash content to queryhash
							$queryhash{$current} = $querycount;
							$queryseqnum[$querycount] = 0;
							foreach $tmp1 (keys %tmphash)
							{
								$queryseq[$querycount][$queryseqnum[$querycount]] = $tmp1;
								$queryseqnum[$querycount]++;
							}
							$querycount++;
						}
					}
					$current = $inlinehmm;
					$querylen[$querycount] = $arr[5];
					%tmphash = ();
				}
				$i = $arr[16] - $arr[15];
				if ($i / $arr[5] >= $COV_CUTOFF)
				{
					if (!(exists $tmphash{$1}))
					{
						$tmphash{$1} = 1;
					}
				}
			}
		}
	}
	if ((scalar(keys %tmphash)) > 0)
	{
		# Flush tmphash content to queryhash
		$queryhash{$current} = $querycount;
		$queryseqnum[$querycount] = 0;
		foreach $tmp1 (keys %tmphash)
		{
			$queryseq[$querycount][$queryseqnum[$querycount]] = $tmp1;
			$queryseqnum[$querycount]++;
		}
		$querycount++;
	}
	close(LOCALFILE);

	if ((scalar keys %queryhash) == 0)
	{
		return(-1);
	}

=print marker gene information
	%tmphash = ();
	foreach $tmp1 (keys %queryhash)
	{
		print "$tmp1 -> $queryseqnum[$queryhash{$tmp1}]\n";
		if (!(exists $tmphash{$queryseqnum[$queryhash{$tmp1}]}))
		{
			$tmphash{$queryseqnum[$queryhash{$tmp1}]} = 1;
		}
		else
		{
			$tmphash{$queryseqnum[$queryhash{$tmp1}]}++;
		}
	}
	foreach $tmp1 (sort {$a <=> $b} keys %tmphash)
	{
		print "$tmp1, $tmphash{$tmp1}\n";
	}
=cut

    @arr = sort {$a <=> $b} (@queryseqnum);
    $i = $arr[int($querycount / 2)];
    print "$i";

	# Check the number of markers
	if ($i <= 1)
	{
        return -1;
	}
	else
	{
		$j = 999999;
		foreach $tmp1 (keys %queryhash)
		{
			if ($queryseqnum[$queryhash{$tmp1}] == $i && $j > $querylen[$queryhash{$tmp1}])
			{
				$j = $querylen[$queryhash{$tmp1}];
				$current = $tmp1;
			}
		}
		#print "$current\n";

		open(LOCALOUT, ">$out");
		for ($i = 0; $i < $queryseqnum[$queryhash{$current}]; $i++)
		{
			print LOCALOUT "$queryseq[$queryhash{$current}][$i]\n";
		}
		close(LOCALOUT);
	}
	return 0;
}


sub checkMarker
{
	if ($_[0] eq "TIGR00388")
	{
		return "TIGR00389";
	}
	elsif ($_[0] eq "TIGR00471")
	{
		return "TIGR00472";
	}
	elsif ($_[0] eq "TIGR00408")
	{
		return "TIGR00409";
	}
	elsif ($_[0] eq "TIGR02386")
	{
		return "TIGR02387";
	}
	else
	{
		return $_[0];
	}
}

sub trim
{
	my $t = $_[0];
	my $i;
	$t =~ s/\s+$//;
	$i = index($t, ' ');
	if ($i >= 0)
	{
		$t = substr($t, 0, $i);
	}
	return $t
}


1;
