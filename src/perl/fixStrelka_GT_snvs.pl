my @nfields;
my @samples;
while(<>)       {
        if (/^\#\#/) {
                print $_;
                next;
        }
        if (/^\#/) {
                print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
                chomp;
                my $line = $_;
                $line =~ s/\#//;
                @nfields = split (/\t/, $line) ;
                @samples = @nfields[9..$#nfields];
                print "#$line\n";
                next;
        }
        chomp;
        my %fields;
        my @values = split/\t/;
        
        foreach my $i (@nfields) {
                        $fields{$i} = shift @values;
        }
        my @names_format = split (/:/, $fields{'FORMAT'});
        push(@names_format, 'GT');
	my @values_info = split (/;/, $fields{'INFO'});
	
	my ($NormalAllele,$TumorAllele);		
	foreach my $i (@values_info)	{
		if ($i =~ /^SGT=/)	{
			my ($sgt,$gtvalue)=split('=',$i);
			($NormalAllele,$TumorAllele)=split('->',$gtvalue);
			$NormalAllele =~ s/$fields{'REF'}/0/g;
			$NormalAllele =~ s/$fields{'ALT'}/1/g;
			$TumorAllele =~ s/$fields{'REF'}/0/g;
                        $TumorAllele =~ s/$fields{'ALT'}/1/g;
			if ($NormalAllele eq "00")	{
				$NormalAllele="0/0";
			}elsif ($NormalAllele eq "01" || $NormalAllele eq "10")	{
				$NormalAllele="0/1";
			}else	{$NormalAllele="1/1";}	
			if ($TumorAllele eq "00")      {
                                $TumorAllele="0/0";    
                        }elsif ($TumorAllele eq "01" || $TumorAllele eq "10") {
                                $TumorAllele="0/1";
                        }else   {$TumorAllele="1/1";}
		}
	}
	my (%sample_values);
	foreach my $sample (@samples) {
		my @values;
		my @s_values = split (/:/, $fields{$sample});
		foreach my $i (@names_format) {
			my $cvalue = shift (@s_values);
			$sample_values{$sample}->{$i}=$cvalue;
			if ($i eq 'GT')	{
				if ($sample eq 'NORMAL')	{
					push(@values, $NormalAllele);
				}else { 
					push(@values,$TumorAllele);
				}
			}else	{
				push(@values, $sample_values{$sample}->{$i});
			}
		}
		@values=reverse @values;
	$fields{$sample} = join (':', @values);
	}
	@names_format= reverse @names_format;	
my @line_fields;
        $fields{'FORMAT'} = join (':', @names_format);
        foreach my $field (@nfields) {
        push (@line_fields, $fields{$field});
    }
    print join ("\t",@line_fields) . "\n";
}





