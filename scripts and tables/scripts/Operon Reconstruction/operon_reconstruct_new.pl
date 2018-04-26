#operon_reconstruct_new.pl - generate a genome proximity network (multi-partite) of putative EPS subunit hits
#Question: which hits/sequence families are frequently found in proximity to each other in bacterial genomes?

#Calculate the average genome proximity (in KBs) between genes encoding putative EPS subunit hits 
#that belong to distinct sequence groups (95% needleman-wunsch pairwise alignment sequence identity).
#Do sequence groups for distinct EPS subunit contain sequences from the same set of genomes?


use strict;
use warnings;

my $infile;
my $outfile;
my $target;

my @liner;

my $sys;
my $subunit;
my %hit_sets;

my $chosen = shift @ARGV;

my @c = ("cellulose", "acet-cellulose", "pel", "pnag", "alginate");

my @ps = ("BcsA", "WssB", "PelF", "PgaC", "Alg8");
my @sub_other = ("Bcs[BZC]", "Wss[FGHI]", "Pel[ABCDEG]", "Pga[ABD]", "Alg[44KEGXLIJFAD]");
my $ps;
my $s_other;

my $q;

for($q = 0; $q <= $#c; $q++){
	if($chosen eq $c[$q]){
		$ps = $ps[$q];
		$s_other = $sub_other[$q];
		last;
	}
}

foreach $infile (<"./eps_hits/$chosen/*_hits_loc.txt">){
	
	if($infile =~ /eps_hits\/(.*)\/(.*)_hits_loc/){
		$sys = $1;
		$subunit = $2;
		if($infile !~ /$sys\/(.*)_(no_cut|1e-2)/){
			$hit_sets{$sys}{$subunit} = $infile;
		}
	}

}

my $group_num;
my %group;
my $gi;


my $sp;
my $sp_id;
my $prot;
my $loc;
my %hits;

my $st;
my $en;


my %sp;

foreach $sys (keys %hit_sets){
	foreach $subunit (keys %{$hit_sets{$sys}}){
		$infile = $hit_sets{$sys}{$subunit};

		open IN, $infile or die;
		while(<IN>){
			chomp;
			@liner = split /\t/;
			($prot, $gi, $sp_id, $sp, $loc) = (@liner[0..1], @liner[2..3], $liner[5]);
			#print join "\t", ($prot, $gi, $sp, $loc), "\n";
			#print $sp_id, "\t", $sp, "\n";

			if($sp_id){
				$sp{$sp_id} = $sp;
			}
			elsif(!$sp_id && $sp){
				$sp{$sp} = $sp_id;
			}
			
			
			#if($sp_id  && $sp_id eq "NC_016027"){
				#print $sp_id, "\t", $sp, "\t", $subunit, "\t", $gi, "\t", $loc, "\n";
			
			if($loc){
				($st, $en) = parse_loc($loc);
				if(!exists $hits{$sys}{$sp_id}{$gi}){
					$hits{$sys}{$sp_id}{$gi} = [$subunit, $gi, $loc] if $sp_id && $loc;
				}
				else{
					$hits{$sys}{$sp_id}{$gi}[0] = join "+", $hits{$sys}{$sp_id}{$gi}[0], $subunit;
				}
			}
			#}
			
			
			#print $_, "\n" if $prot =~ /YP_008047743|YP_112018/;
		}
		close IN;
	}
}

my $t1;
my $t2;
my @t1;
my @t2;
my $i;
my $j;
my $k;

my @g1;
my @g2;
my @p1;
my @p2;

my $g1;
my $g2;
my $p1;
my $p2;

my $sp1;
my $sp2;
my $l1;
my $l2;

my $s1;
my $e1;
my $s2;
my $e2;

my $dist;

my %hit_dist;
my %group_tally;

my @sys;
my $sys1;
my $sys2;

#for each pair of targets / eps subunits (t1, t2)
#for each pair of target groups (g1, g2)
#for each pair of proteins belonging to the same species (p1, p2)

my %crap_stack;

my @sub;
my $sub1;
my $sub2;
my @hits;

my $prot1;
my $prot2;

my %op;

my $op_num;

my %dist_cuts = ("acet-cellulose" => 7500,
				"cellulose" => 4000,
				"pel" => 10000,
				"pnag" => 2000,
				"alginate" => 7500
);

my $dist_cut = $dist_cuts{$chosen};

my %chosen;
my $next = 0;

foreach $sys (keys %hits){
	foreach $sp_id (keys %{$hits{$sys}}){
	#if($sp_id eq "NC_016935"){
		@hits = sort {$hits{$sys}{$sp_id}{$a}[2] cmp $hits{$sys}{$sp_id}{$b}[2]} keys %{$hits{$sys}{$sp_id}};
		$op_num = 1;
		#print join "\n", @hits, "\n";
		for($i = 0; $i < $#hits; $i++){
			($sub1, $prot1, $l1) = @{$hits{$sys}{$sp_id}{$hits[$i]}};
			
			($s1, $e1) = parse_loc($l1);
			$crap_stack{$hits[$i]} = join "\t", $sub1, $sp1, $p1, $l1 if !$e1;
		
			
				
				($sub2, $prot2, $l2) = @{$hits{$sys}{$sp_id}{$hits[$i+1]}};
				($s2, $e2) = parse_loc($l2);
				
				
				
				if($sys ne "acet-cellulose"){
					$chosen{$sp_id} = 1 if $sub1 =~ /$ps/ || $sub2 =~ /$ps/;
				}
				else{
					$chosen{$sp_id} = 1 if $sub1 =~ /$ps/ || $sub2 =~ /$ps/;
				}
				#$chosen{$sp_id} = 1 if $sub1 =~ /Wss(FGHI)/ || $sub2 =~ /Wss(FGHI)/;
				
				$crap_stack{$hits[$j]} = join "\t", $sub2, $sp2, $p2, $l2 if !$e2;
				
				#print $l1, "\t", $l2, "\n";

				$dist = $s1 - $e2 if $s1 > $e2;
				$dist = $s2 - $e1 if $s1 < $e2;
				
				print $sub1, "\t", $sub2, "\t", $l1, "\t", $l2, "\t", $dist, "\n" if $sp_id eq "NC_009662";
				
				if($l1 =~ /complement/ && $l2 !~ /complement/){
				#	$op_num++;
				#	next;
				}
				elsif($l1 !~ /complement/ && $l2 =~ /complement/){
				#	$op_num++;
				#	next;
				}
				
				if($sp_id eq "NC_009952" || $sp_id eq "NC_009662"){
					$dist_cut = 10000;
				}
				else{
					$dist_cut = $dist_cuts{$chosen};
				}
				
				if($dist <= $dist_cut){
					$op{$sys}{$sp_id}{$op_num}{$l1} = join "|", $sub1, $prot1, $l1;
					$op{$sys}{$sp_id}{$op_num}{$l2} = join "|", $sub2, $prot2, $l2;
				}
				else{
					$op_num++;
				}
		}
	#}
	}
}



my $avg_dist;
my $num_sp;
my $group1;
my $group2;

my $sp_name;
my $tax;
my $op_size;

$outfile = "./$chosen\_op_final.txt";

print $outfile, "\n";

open OUT, ">", $outfile;

my $op;

#print out reconstructed operons!!!!

foreach $sys (keys %op){

	#foreach $sp_id (sort {keys %{$op{$sys}{$b}} <=> keys %{$op{$sys}{$a}}} keys %{$op{$sys}}){
	foreach $sp_id (sort {$a cmp $b} keys %chosen){
		#if($sp_id eq "NC_009952"){
		foreach $op_num (sort {$a cmp $b} keys %{$op{$sys}{$sp_id}}){
			$sp = $sp{$sp_id};

			
			($sp_name, $tax) = split /\| /, $sp;
			#if($tax =~ /bacillus/i){
			$op_size = keys %{$op{$sys}{$sp_id}{$op_num}};
		
			if($op_size > 1){
				$op = "";
				
		
				foreach $loc (sort {$a cmp $b} keys %{$op{$sys}{$sp_id}{$op_num}}){
					($subunit, $gi) = split /\|/, $op{$sys}{$sp_id}{$op_num}{$loc};
					$prot = $hits{$sys}{$subunit}{$gi}[2];
					#print "\t", join "|", $subunit, $prot, $loc;
					$op = join "\t", $op, (join "|", $subunit, $gi, $loc);
					
	
				}
				
				if($op =~ /$ps/ && $op =~ /$s_other/){
					#print $sp_name, "\t", $tax, "\t", $sp_id, "\t", $op_size, $op, "\n";
					print OUT $sp_name, "\t", $tax, "\t", $sp_id, "\t", $op_size, $op, "\n";
				}
				elsif($op =~ /WssB/ && $op =~ /Wss[FGHI]/){
					#print OUT $sp_name, "\t", $tax, "\t", $sp_id, "\t", $op_size, $op, "\n";
				}
				#print OUT "\t", join "|", $subunit, $gi, $loc;

			}
		}
		#}
	}
}

close OUT;

#Check for subunit hit pairs for which a genomic distance could not be calculated
foreach $p1 (keys %crap_stack){
	#print $crap_stack{$p1}, "\n";

}

sub parse_loc{
	my $loc = shift @_;
	my $l = $loc;
	my $s;
	my $e;
	
	$l =~ s/complement//;
	$l =~ s/\(|\)//g;
	$l =~ s/\<|\>//;
	
	$l =~ s/join//;
	
	if($loc =~ /complement/){
		($s, $e) = split /\.+/, $l if $loc !~ /join/;
		
		if($loc =~ /join/){
			($s, $e) = split /,/, $l;
			$s =~ s/.*\.+//;
			$e =~ s/\.+.*//;
		}
	}
	else{
		($s, $e) = split /\.+/, $l if $loc !~ /join/;
		
		if($loc =~ /join/){
			($s, $e) = split /,/, $l;
			$s =~ s/\.+.*//;
			$e =~ s/.*\.+//;
		}
		
	}
	
	return ($s, $e);
	
}
