package utils;
#*************************************************************************
#
#   Program:    
#   File:       utils.pm
#   
#   Version:    V1.4
#   Date:       26.06.19
#   Function:   
#   
#   Copyright:  (c) Prof. Andrew C. R. Martin, UCL, 2015-2019
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   General utility functions
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   01.05.15 Original   By: ACRM
#   V1.1   17.05.16 Added intellisplit()
#   V1.2   26.06.19 Added various functions from abYmod util.pm
#   V1.3   15.08.19 Added CreateTempFileName()
#   V1.4   25.07.21 BuildPackage now takes .tjz and .tar.gz as well
#                   as .tgz
#
#*************************************************************************
use File::Basename;

#*************************************************************************
# Prints a message ($string) with an optional line number ($line) and 
# exits the program
sub mydie
{
    my($string, $line) = @_;
    print STDERR "$string";
    if($line)
    {
        print STDERR " at line $line";
    }
    print "\n";
    exit 1;
}

#-------------------------------------------------------------------------
# Installs $file in $dir and uncompresses it if $uncompress is set
sub InstallFile
{
    my($file, $dir, $uncompress) = @_;

    return(0) if(!MakeDir($dir));
    `cp $file $dir`;

    my $newfile = $file;
    $newfile =~ s/.*\///;          # Remove path
    $newfile = "$dir/$newfile";
    return(0) if(! -e $newfile);

    if($uncompress)
    {
        `cd $dir; gunzip -f $newfile`;
        $newfile =~ s/\.gz//;
        return(0) if(! -e $newfile);
    }

    return(1);
}

#-------------------------------------------------------------------------
# Makes a directory and checks that it has been created OK
sub MakeDir
{
    my($dir) = @_;
    `mkdir -p $dir` if(! -d $dir);
    return(0) if(! -d $dir);
    return(1);
}

#-------------------------------------------------------------------------
# Checks that a list of files are executable. Returns a list of the files
# that were NOT OK, or a blank string if all were OK. The list may also
# be supplied as a scalar variable, separating names with |
sub CheckExecutables
{
    my(@files) = @_;

    # If there is only one filename specified and it contains a |
    # then split this into a list
    if((scalar(@files) == 1) && ($files[0] =~ '\|'))
    {
        @files = split(/\|/, $files[0])
    }

    my $badexe = '';
    foreach my $file (@files)
    {
        if(! -x $file)
        {
            if($badexe eq '')
            {
                $badexe = $file;
            }
            else
            {
                $badexe .= ", $file";
            }
        }
    }
    return($badexe);
}

#-------------------------------------------------------------------------
# Checks that a list of environment variables have been defined. Returns 
# a list of those that were NOT OK, or a blank string if all were OK. The 
# list may also be supplied as a scalar variable, separating names with |
sub CheckEnvironmentVariables
{
    my(@vars) = @_;

    # If there is only one filename specified and it contains a |
    # then split this into a list
    if((scalar(@vars) == 1) && ($vars[0] =~ '\|'))
    {
        @vars = split(/\|/, $vars[0])
    }

    my $badvar = '';
    foreach my $var (@vars)
    {
        if(!defined($ENV{$var}))
        {
            if($badvar eq '')
            {
                $badvar = $var;
            }
            else
            {
                $badvar .= ", $var";
            }
        }
    }
    return($badvar);
}

#-------------------------------------------------------------------------
# Checks that a list of files exist and are readable. Returns a list of 
# the files that were NOT OK, or a blank string if all were OK.
sub CheckFile
{
    my(@files) = @_;
    my $badfile = '';
    foreach my $file (@files)
    {
        if(! -r $file)
        {
            if($badfile eq '')
            {
                $badfile = $file;
            }
            else
            {
                $badfile .= ", $file";
            }
        }
    }
    return($badfile);
}


#-------------------------------------------------------------------------
sub parseFilename
{
    my($infile, $longext) = @_;
    my ($path, $filename, $filestem, $extension) = ('','','','');

#OK
    $filename = $infile;
    $filename =~ s/^.*\///;     # Remove anything up to the first /

#OK
    if($longext)
    {
        if($filename =~ /\.(.+)$/)
        {
            $extension = $1;
        }
    }
    else
    {
        if($filename =~ /.*\.(.+?)$/)
        {
            $extension = $1;
        }
    }

# OK
    if($infile =~ /(.*)\//)
    {
        $path = $1;
    }

    if($longext)
    {
        if($filename =~ /^(.*?)\..*/)
        {
            $filestem = $1;
        }
    }
    else
    {
        if($filename =~ /^(.*)\./)
        {
            $filestem = $1;
        }
    }

    return($path, $filename, $filestem, $extension);
}

#-------------------------------------------------------------------------
sub setExtension
{
    my($infile, $ext, $longext) = @_;
    my $outfile;

    if(!($ext =~ /^\./))        # If it doesn't start with a . add one
    {
        $ext = ".$ext";
    }

    my($path, $filename, $filestem, $extension) = parseFilename($infile, $longext);
    if($extension eq '')
    {
        $outfile = $infile . $ext;
    }
    else
    {
        if($path eq '')
        {
            $outfile = "$filestem$ext";
        }
        else
        {
            $outfile = "$path/$filestem$ext";
        }
    }
    return($outfile);
}

#-------------------------------------------------------------------------
sub ReadFileHandleAsTwoColumnHashRef
{
    my ($fh)   = @_;
    my %result = ();

    while(my $line = <$fh>)
    {
        chomp $line;
        $line =~ s/\#.*//;  # Remove comments
        $line =~ s/^\s+//;  # Remove leading spaces
        $line =~ s/\s+$//;  # Remove trailing spaces
        if(length($line))
        {
            my @fields = split(/\s+/, $line);
            $result{$fields[0]} = $fields[1];
        }
    }

    return(\%result);
}

#-------------------------------------------------------------------------
# @fields = intellisplit($string)
# -------------------------------
# Like split but allows single or double inverted commas to wrap a string
# including spaces. (Double inverted commas may be contained in a pair
# of single inverted commas and vice versa.) It only splits at a normal 
# space, not at tabs.
#
# 17.05.16 Original   By: ACRM

sub intellisplit
{
    my($input) = @_;
    my @in = split(//, $input);
    my $inDic = 0;
    my $inSic = 0;
    my @output = ();
    my $string = '';

    foreach my $char (@in)
    {
        if($char eq '"')
        {
            if(!$inSic)
            {
                $inDic = $inDic?0:1;
            }
            $string .= $char;
        }
        elsif($char eq "'")
        {
            if(!$inDic)
            {
                $inSic = $inSic?0:1;
            }
            $string .= $char;
        }
        elsif(($char eq ' ') && !($inDic || $inSic))
        {
            push @output, $string;
            $string = '';
        }
        else
        {
            $string .= $char;
        }
    }
    push @output, $string;

    return(@output);
}

#*************************************************************************
# my ($filename, $stem) = utils::BuildFileName($inFile, $dir, $ext)
# -----------------------------------------------------------------
# Input:   text  $inFile    Input filename (maybe with a path)
#          text  $dir       Required path
#          text  $ext       Required extension
# Returns  text  $filename  New filename
#          text  $stem      The filestem
#
# Builds a filename by discarding the path and extension from an input 
# filename and adding new path and extension
#
# For some odd reason /\..*?$/ doesn't properly do a non-greedy match -
# it matches from the first '.' instead of the last '.'. Consequently 
# the code to remove te extension has to do a greedy match on the first
# part of the string and substitute that for the whole thing
#
#  19.09.13  Original  By: ACRM
sub BuildFileName
{
    my($inFile, $dir, $ext) = @_;

    my $stem = $inFile;
    $stem =~ s/(.*)\..*?$/$1/;           # Remove extension
    $stem =~ s/.*\///;                   # Remove path

    chop $dir if($dir =~ /\/$/);         # Remove / from end of path
    $ext = ".$ext" if(!($ext =~ /^\./)); # Prepend . to extension if needed

    my $outFile = "$dir/$stem$ext";      # Construct filename

    return($outFile, $stem);
}

#*************************************************************************
#> $tmpDir = CreateTempDir($pName)
#  -------------------------------
#  Input:   string   $pName    Directory base name
#  Return:  string             Full created directory name
#
#  Create a temporary directory. The time is appended to the filestem
#  and this is placed in the directory named in $config::tmp
#
#  19.09.13  Original  By: ACRM
sub CreateTempDir
{
    my ($pName) = @_;
    my $tmpDir = $config::tmp . "/${pName}_$$" . "_" . time;
    $tmpDir =~ s/\/\//\//g; # Replace // with /
    `mkdir -p $tmpDir`;
    if(! -d $tmpDir)
    {
        return(undef);
    }
    return($tmpDir);
}

#*************************************************************************
# @files = GetFileList($dir, $type, $prepend)
# -------------------------------------------
# Input:   string   $dir      A directory path
#          string   $type     A file extension (or blank)
#          BOOL     $prepend  Prepend the directory path onto the filenames
# Returns: string[]           List of filenames
#
# This function gets a list of filenames from a directory which have the
# specified extension - strictly this is just text that the filename must
# contain - if it really must be an extension then it needs to end with a $
# (end of string marker).
# By default, the filenames are returned without the path information.
# If $prepend is set, then the path is prepended onto each filename
#
#  19.09.13  Original  By: ACRM
sub GetFileList
{
    my($dir, $type, $prepend) = @_;
    my @files = ();

    chop $dir if($dir =~ /\/$/);

    if(opendir(my $dh, $dir))
    {
        @files = grep(!/^\./, readdir($dh));
        if($type ne "")
        {
            @files = grep(/$type/, @files);
        }
        closedir($dh);
    }
    if($prepend)
    {
        foreach my $file (@files)
        {
            $file = "$dir/$file";
        }
    }
    return(sort(@files));
}

#*************************************************************************
#> BOOL FileNewer($testFile, $refFile)
#  -----------------------------------
#  Input:   string   $testFile   File of which to check date
#           string   $refFile    Reference file against which we compare
#  Returns: BOOL                 True if either file does not exist or
#                                if $testFile is newer than $refFile
#
#  Tests whether $testFile is newer than $refFile
#
#  19.09.13  Original  By: ACRM
sub FileNewer
{
    my($testFile, $refFile) = @_;

    return(1) if((! -e $testFile) || (! -e $refFile));

    my @stats;
    @stats = stat($testFile);
    my $testDate = $stats[9];
    @stats = stat($refFile);
    my $refDate = $stats[9];

    return((($testDate > $refDate)?1:0));
}

#*************************************************************************
#> void CheckAndDie($path, $isDir, $text)
#  --------------------------------------
#  Input:   string   $path   Path to a directory or file
#           BOOL     $isDir  $path is a directory
#           string   $text   Text message
#
#  Tests whether $path is a directory or file (depending on $isDir) and
#  if it isn't what it is supposed to be, prints a message including $text
#  and dies
#
#  19.09.13  Original  By: ACRM
sub CheckAndDie
{
    my($path, $isDir, $text) = @_;
    my $ok = 1;
    if($isDir)
    {
        $ok = 0 if(! -d $path);
    }
    else
    {
        $ok = 0 if(! -f $path);
    }

    if(!$ok)
    {
        mydie("\nabYmod configuration/installation error: $text doesn't exist:\n   $path\n\n");
    }
}

#*************************************************************************
#> %hash = BuildTwoColumnHash(@array)
#  ----------------------------------
#  Input:   string[]   @array    An array of strings of the form 'X Y'
#  Returns: hash                 A hash where X is the key and Y the value
#
#  Builds a hash from an array containing strings with two columns. 
#  NOTE! - If items in the first column (X) are repeated, then the last occurrence
#          will be used as the stored value
#        - If the string contains more than two columns, then only the second
#          column will be stored (i.e. 'X Y Z' will use X as a key for Y and Z
#          will be discarded.
#
#
#  19.09.13  Original  By: ACRM
sub BuildTwoColumnHash
{
    my @records = @_;
    my %resultHash = ();
    foreach my $line (@records)
    {
        chomp $line;
        my @fields = split(/\s+/, $line);
        $resultHash{$fields[0]} = $fields[1];
    }
    return(%resultHash);
}

#*************************************************************************
#> BOOL inlist($item, @list)
#  -------------------------
#  Input:   string   $item    An item for which to search
#           string[] @list    An array
#  Return:  BOOL              Was $item in the array?
#
#  Tests whether an item appears in the array
#
#  19.09.13  Original  By: ACRM
sub inlist
{
    my($item, @list) = @_;
    foreach my $listItem (@list)
    {
        if($item eq $listItem)
        {
            return(1);
        }
    }
    return(0);
}

#*************************************************************************
# @textArray = ReadFileAsArray($inFile)
# -------------------------------------
# Input:   text   $inFile    A file name to be read
# Returns: text[]            An array of lines from the file
#
#  19.09.13  Original  By: ACRM
sub ReadFileAsArray
{
    my($inFile) = @_;
    my @contents = ();

    if(open(my $fp, $inFile))
    {
        while(<$fp>)
        {
            chomp;
            s/\r//;
            push @contents, $_;
        }
        close $fp;
    }
    return(@contents);
}

#*************************************************************************
#> void RunCommand($exe, $useSystem)
#  ---------------------------------
#  Input:   string  $exe    An excutable string
#
#  Runs a command
#  19.09.13  Original  By: ACRM
#  28.09.15  Now returns the output
#  30.09.15  Added $useSystem parameter (optional)
sub RunCommand
{
    my ($exe, $useSystem) = @_;
    my $result = '';

    print STDERR "$exe\n";
    if(defined($useSystem) && $useSystem)
    {
        system("$exe");
    }
    else
    {
        $result = `$exe`;
    }
    return($result);
}

#*************************************************************************
#> $dir = GetDir($file)
#  --------------------
#  Input:   string   $file    Full path to a file
#  Return:  string            The directory (path) element of the filename
#
#  19.09.13  Original  By: ACRM
sub GetDir
{
    my ($file) = @_;
    $file =~ /(.*\/)/;  
    my $dir = $1;
    return($dir);
}

#*************************************************************************
#> @result = sortArrayByArray($aTarget, $aKey)
#  -------------------------------------------
#  Input:   \data[]  $aTarget  Reference to array to be sorted
#           \data[]  $aKey     Reference to array on which to sort
#  Returns: data[]   @result   Sorted version of $aTarget
#
#  Sorts the $aTarget array based on the values in the $aKey array
#
#  17.07.14 Original   By: ACRM
sub sortArrayByArray
{
    my ($aTarget, $aKey) = @_;
    my @idx = sort {$$aKey[$a] <=> $$aKey[$b]} 0 .. $#$aKey;
    my @target = @$aTarget[@idx];
    return(@target);
}

#*************************************************************************
#> void DEBUG($string)
#  -------------------
#  Input:   string   $string   A text string
#
#  Prints the string if $::debug is defined
#
#  17.07.14  Original   By: ACRM
sub DEBUG
{
    my($string) = @_;

    if(defined($::debug))
    {
        print STDERR "DEBUG: $string\n";
    }
}

#*************************************************************************
#>void PrettyPrint($fp, $sequence, $width, $append)
# -------------------------------------------------
# Input:   FILE   $fp         Reference to file handle
#          string $sequence   Sequence to be printed
#          int    $width      Width to print
#          string $append     String to be appended to the sequence
#
# Prints a string breaking it up into $width chunks. The $append string
# is appended to the main string before this is done.
#
# 21.07.14 Original   By: ACRM
sub PrettyPrint
{
    my($fp, $sequence, $width, $append) = @_;
    $sequence .= $append;

    while($sequence ne '')
    {
        print $fp substr($sequence, 0, $width) . "\n";
        $sequence = substr($sequence, $width);
    }
}


#*************************************************************************
#> void BuildPackage($package, $subdir, $aExe, $binDir, $dataDir, 
#                    $dataDest, $postBuild)
#  -------------------------------------------------------------------------
#  Input:   string  $package    The gzipped tar file of the package
#           string  $subdir     Subdirectory of the unpacked package
#                               containing source code
#           string  $aAxe       Reference to array of exectuables generated
#           string  $binDir     Destination binary directory
#           string  $dataDir    Data directory in unpacked package
#           string  $dataDest   Destination data directory
#
#  Builds and installs a C package
#
#  19.09.13  Original  By: ACRM
#  25.09.15  Now takes a reference to an array of executables
#  01.10.15  Makes the destination directories if they don't exist
#  02.10.15  Moved into util.pm
#  13.09.16  Added checks that install of executable files has 
#            actually worked
#  24.03.17  Added $postBuild
sub BuildPackage
{
    my ($package, $subdir, $aExe, $binDir, $dataDir, $dataDest, $postBuild) = @_;

    # See if we need to do this - i.e. we don't have the files
    # already
    my $needsToRun = 0;
    foreach my $exe (@$$aExe)
    {
        if(! -x "$binDir/$exe")
        {
            $needsToRun = 1;
            last;
        }
    }
    if(($dataDir ne '') && ( ! -d $dataDest))
    {
        $needsToRun = 1;
    }

    if($needsToRun)
    {
        utils::RunCommand("tar xvf $package");
        my $packageDir = $package;
        $packageDir =~ s/.*\///;
        $packageDir =~ s/\.tgz//;
        $packageDir =~ s/\.tar\.gz//;
        $packageDir =~ s/\.tjz//;
        utils::RunCommand("cd $packageDir/$subdir; make");
        foreach my $exe (@$$aExe)
        {
            utils::RunCommand("mkdir -p $binDir") if(! -d $binDir);
            utils::RunCommand("cp $packageDir/$subdir/$exe $binDir");

            if(! -e "$binDir/$exe")
            {
                mydie("\nabYmod installation error: $binDir/$exe not created.\n       Compilation in $packageDir probably failed\n\n");
            }
        }
        if($dataDir ne "")
        {
            utils::RunCommand("mkdir -p $dataDest") if(! -d $dataDest);
            utils::RunCommand("cp -R $packageDir/$dataDir/* $dataDest");
        }
        if($postBuild ne "")
        {
            utils::RunCommand($postBuild);
        }
        `rm -rf $packageDir`;
    }
    else
    {
        print STDERR "Skipped installation of $package - already installed\n";
    }
}


#*************************************************************************
#> BOOL CheckLibrary(@libs)
#  ------------------------
#  Input:    @libs   Array of library names to search for
#  Returns:  BOOL    Found?
#
#  Checks the library paths to see if a specified library exists
#
#  13.09.16 Original   By: ACRM
sub CheckLibrary
{
    my(@files) = @_;

    my @dirs = qw(/usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64);

    return(CheckFilesExistInDirs(\@files, \@dirs));
}


#*************************************************************************
#> BOOL CheckInclude(@incs)
#  ------------------------
#  Input:    @libs   Array of include file names to search for
#  Returns:  BOOL    Found?
#
#  Checks the system include paths to see if a specified library exists
#
#  13.09.16 Original   By: ACRM
sub CheckInclude
{
    my(@files) = @_;

    my @dirs = qw(/usr/include /usr/local/include);

    return(CheckFilesExistInDirs(\@files, \@dirs));
}


#*************************************************************************
# BOOL CheckFilesExistInDirs(\@FilesToCheck, \@DirsToSearch)
# ----------------------------------------------------------
# Input:    \@FilesToCheck  Ref to list of files we are looking for
#           \@DirsToSearch  Ref to list of directories to search
#
# Checks if the specified files exist in any of the specified directories.
# These are searched recursively.
#
#  13.09.16 Original   By: ACRM
sub CheckFilesExistInDirs
{
    my($aFiles, $aDirs) = @_;

    foreach my $inFile (@$aFiles)
    {
        my $found = 0;

        foreach my $location (@$aDirs)
        {
            if(!$found)
            {
                my @fileList = GetRecursiveFileList($location);

                foreach my $file (@fileList)
                {
                    
                    if($file =~ /$inFile$/)
                    {
                        $found = 1;
                        last;
                    }
                }
            }
        }

        if(!$found)
        {
            return(0);
        }
    }

    return(1);

}


#*************************************************************************
#> @files = GetRecursiveFileList($location)
#  ----------------------------------------
#  Input:    $location   Top level directory
#  Returns:              List of full file paths in that directory
#
#  Builds a list of all files below a given direcotory - uses 'ls -R' to
#  obtain a recursive list
#
#  13.09.16 Original   By: ACRM
sub GetRecursiveFileList
{
    my($location) = @_;
    my $stem = '';
    my @files = ();
    my $dirTree = `\\ls -R $location 2>/dev/null`;
    my @records = split(/\n/, $dirTree);
    foreach my $record (@records)
    {
        $record =~ s/\s+$//;    # Remove trailing whitespace
        if(length($record))
        {
            if($record =~ /(.*)\:/)
            {
                $stem = $1 . '/';
            }
            else
            {
                push @files, "$stem$record";
            }
        }
    }

    return(@files);
}

#*************************************************************************
#> ($bad, $msg) = CheckModulesInstalled(@modules)
#  ----------------------------------------------
#  Input:    @modules   Array of module names
#  Returns:  $bad       >0: the number of missing modules
#                        0: all modules found
#            $msg       Names of missing modules
#
#  Checks whether given mmodules are installed.
#
#  28.06.19 Original   By: ACRM
sub CheckModulesInstalled
{
    my(@modules) = @_;
    
    my $msg = '';
    my $bad = 0;
    for(@modules)
    {
        eval "use $_";
        if($@)
        {
            $msg .= "-  $_\n";
            $bad++;
        }
    }

    return($bad, $msg);
}

#*************************************************************************
sub CreateTempFileName
{
    my ($stem) = @_;
    my $fnm = "/var/tmp/${stem}_" . $$ . time();
    return($fnm);
}

        
1;
