my @nfields;
my @samples;
while(<>)	{
	if (/^\#\#/) {
		print $_;
		next;
	}
	if (/^\#/) {
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
 	push(@names_format, 'DP');
	my @values_info = split (/;/, $fields{'INFO'});
	my (%sample_values);
	foreach my $sample (@samples) {
		my @values;
		my @s_values = split (/:/, $fields{$sample});
		foreach my $i (@names_format) {
			my $cvalue = shift (@s_values);
			$sample_values{$sample}->{$i}=$cvalue;
			if ($i eq 'DP')	{
				if ($sample_values{$sample}-> {'GT'} ne './.')	{
					($refV,$refAlt)=split(/,/,$sample_values{$sample}-> {'AD'});
					$DPT=$refV+$refAlt;
					$sample_values{$sample}-> {'DP'}=$DPT;
					push(@values, $sample_values{$sample}-> {'DP'});
				}else	{
					push(@values, $sample_values{$sample}->{$i});
				}
			} else	{
				if ($sample_values{$sample}-> {'GT'} eq './.')	{
					@values=('./.');
				} else { push(@values, $sample_values{$sample}->{$i});}
			}
			
		}
		 $fields{$sample} = join (':', @values);
			if ($fields{$sample} eq './.:')	{
				$fields{$sample}="./.";
			}	
	}
	my @line_fields;
    	$fields{'FORMAT'} = join (':', @names_format);
	foreach my $field (@nfields) {
        push (@line_fields, $fields{$field});
    }
    print join ("\t",@line_fields) . "\n";
}

	
