#! /usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage; 
use List::Util 'max';
use List::Util 'min';
my $flanklength=35;
my $grnafile='';
my $tctsfasta='';
my $message_text  = "please specify option parameters. -i < fasta file containing a set of gRNAs > -t <fasta file containing genes being targeted>";
  my $exit_status   = 2;          ## The exit status to use
  my $verbose_level = 0;          ## The verbose level to use
  my $filehandle    = \*STDERR;   ## The filehandle to write to
GetOptions ('i=s' => \$grnafile, 't=s' => \$tctsfasta);
pod2usage( { -message => $message_text ,
               -exitval => $exit_status  ,  
               -verbose => $verbose_level,  
               -output  => $filehandle } ) if ($grnafile eq '' or $tctsfasta eq '');

my $dirname = "flankseq";
if (-e "$dirname") {
}
else {
	mkdir("$dirname",0777) or die;
}			   

#######read in all TcTS, store all TcTS for finding gRNA targets,
my %allTcTS;
my $seqio=Bio::SeqIO->new(-format => 'Fasta', -file=>"$tctsfasta");
while(my $seqobj=$seqio->next_seq)
{
	my $id=$seqobj->id;
	my $seq=$seqobj->seq;
	$allTcTS{$id}=$seq;
	#print "$seq\n"
}
		   		   
my $seqiog=Bio::SeqIO->new(-format => 'Fasta', -file=>"$grnafile");

while(my $seqobj=$seqiog->next_seq)
{
	my $id=$seqobj->id;
    my $gRNAseq=$seqobj->seq;
	
	#generate result file
	open OUTL, ">flankseq/${id}.Lflank.fasta";
	open OUTR, ">flankseq/${id}.Rflank.fasta";
	$gRNAseq=~s/\|//;
	my $revcom_seqobj = revcom( $gRNAseq );
	my $gRNAseq_revcom=$revcom_seqobj->seq();
	##find on target in new defined TcTS collection
	foreach my $keys (sort keys %allTcTS)
	{
		#print $gRNAseq."\n";
		if ($allTcTS{$keys}=~/$gRNAseq/i) # gRNA is present 
		{
			print OUTL ">".$keys."\n";
			print OUTR ">".$keys."\n";
			for (my $i=0; $i<=$#-; $i++)
			{
				my $Lstart= max (0,$-[$i]-$flanklength+18) ; 
				my $Lend= $-[$i]+18; # leftflank end is the left side of gRNA +18
				my $Lflank= substr ($allTcTS{$keys},$Lstart,$Lend-$Lstart);
				if (length($Lflank)<$flanklength) {$Lflank=('-'x($flanklength-length($Lflank))).$Lflank}
				print OUTL "$Lflank\n";
				
				my $Rstart= $+[$i]-9; # start is the right side of gRNA-9 
				my $Rend= (min (length($allTcTS{$keys}),$+[$i]+$flanklength))-9;
				my $Rflank= substr ($allTcTS{$keys},$Rstart,$Rend-$Rstart);
				if (length($Rflank)<$flanklength) {$Rflank=$Rflank.('-'x($flanklength-length($Rflank)))}
				print OUTR "$Rflank\n";
			}				 
		}
		if ($allTcTS{$keys}=~/$gRNAseq_revcom/i ) # gRNA is present revcom strand
	    {
			#revcom TcTS
			my $revcom_seqobj = revcom( $allTcTS{$keys} );
            my $TcTSseq_revcom=$revcom_seqobj->seq();
			if ($TcTSseq_revcom=~/$gRNAseq/i) # gRNA is present 
			{
				print OUTL ">".$keys."\n";
				print OUTR ">".$keys."\n";
				#print OUT $
				for (my $i=0; $i<=$#-; $i++)
				{
				my $Lstart= max (0,$-[$i]-$flanklength+18) ; 
				my $Lend= $-[$i]+18; # end is the left side of gRNA +18
				my $Lflank= substr ($TcTSseq_revcom,$Lstart,$Lend-$Lstart);
				if (length($Lflank)<$flanklength) {$Lflank=('-'x($flanklength-length($Lflank))).$Lflank}
				print OUTL "$Lflank\n";
				
				my $Rstart= $+[$i]-9; # start is the right side of gRNA -9
				my $Rend= (min (length($allTcTS{$keys}),$+[$i]+$flanklength))-9;
				my $Rflank= substr ($TcTSseq_revcom,$Rstart,$Rend-$Rstart);
				if (length($Rflank)<$flanklength) {$Rflank=$Rflank.('-'x($flanklength-length($Rflank)))}
				print OUTR "$Rflank\n";
				}					 
			} 	 					 
		}						
	}			
	close OUTL;
	close OUTR;	
}