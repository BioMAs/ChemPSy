#!/softs/ActiveTCL/8.6-beta-2/bin/wish8.6
set CTDGenes2DiseasesFile     [lindex $argv 0]
set CTDChemicals2DiseasesFile [lindex $argv 1]
set CTDDiseaseOBOFile         [lindex $argv 2]
set CTDChemicalsFile          [lindex $argv 3]
set EDDir                     [lindex $argv 4]
set HomoloGeneFile            [lindex $argv 5]
set OUTDir                    [lindex $argv 6]

proc HomoloGene4R {HomoloGeneFile OUTDir} {
    #######################################################
    set OUTFile "$OUTDir/HomoloGene4R.txt"
    ::File::Open $OUTFile 100 "no"
    puts "$HomoloGeneFile is loading..."
    set F [open $HomoloGeneFile]
    while {[gets $F Line]>=0} {
	set GeneID ""
	foreach {HGID TAXID GeneID GeneSymb} [lrange [split $Line "\t"] 0 3] {}
	if {$GeneID == ""} {continue}
	set lNewLine [list $GeneID $GeneSymb $TAXID $HGID]
	::File::Save $OUTFile [join $lNewLine "\t"]
    }
    close $F
    puts "$HomoloGeneFile loaded!"
    ::File::Close $OUTFile
    #######################################################

    return $OUTFile
}
proc CTD_chemicals_diseases {CTDChemicalsDiseasesTSVFile OBOCTDDiseasesFile CTDChemicalsFile EDDir OUTDir} {
    set lDiseaseIDs  {}
    set lChemicalIDs {}

    #######################################################
    puts "$CTDChemicalsFile is loading..."
    set F [open $CTDChemicalsFile]
    while {[gets $F Line]>=0} {
	if {[regexp {^\#} $Line]} {continue}
	set CasRNs ""
	foreach {ChemicalName ChemicalID CasRNs} [lrange [split $Line "\t"] 0 2] {}
	if {$CasRNs == ""} {continue}
	regsub {^MESH\:} $ChemicalID "" ChemicalID

	foreach CasRN [split $CasRNs "|"] {
	    if {![info exists TabCas2ChemicalIDs($CasRN)]} {set TabCas2ChemicalIDs($CasRN) {}}
	    lappend TabCas2ChemicalIDs($CasRN) $ChemicalID
	}
    }
    close $F
    puts "$CTDChemicalsFile is loaded!"
    #######################################################

    #######################################################
    foreach EDCASFile [glob -nocomplain "$EDDir/ED_*_CAS.txt"] {
	set DiseaseID [file tail $EDCASFile]
	regsub {_CAS\.txt$} $DiseaseID "" DiseaseID
	lappend lDiseaseIDs $DiseaseID

	puts "$EDCASFile is loading..."
	set F [open $EDCASFile]
	while {[gets $F Line]>=0} {
	    set CasRN [lindex [split $Line "\t"] 0]
	    if {![info exists TabCas2ChemicalIDs($CasRN)]} {continue}
	    foreach ChemicalID $TabCas2ChemicalIDs($CasRN) {
		lappend lChemicalIDs $ChemicalID
		set Tab($ChemicalID,$DiseaseID) 1
	    }
	}
	close $F
	puts "$EDCASFile is loaded!"
    }
    #######################################################

    #########################################
    set lDs {}
    ::OBO::Load   $OBOCTDDiseasesFile
    foreach INIOBOID [ListOfNonRedundantElement [::OBO::Query $OBOCTDDiseasesFile Term]] {
	set TabOBOID2OBOIDWithParents($INIOBOID) {}
	foreach OBOID [ListOfNonRedundantElement [::OBO::Parents $OBOCTDDiseasesFile $INIOBOID 1]] {
	    set OBOName [lindex [OBO::Query $OBOCTDDiseasesFile Term $OBOID name] 0]
	    lappend TabOBOID2OBOIDWithParents($INIOBOID) "$OBOName (${OBOID})"
	    lappend lDs "$OBOName (${OBOID})"
	}
    }
    ::OBO::Delete $OBOCTDDiseasesFile
    set lDiseaseIDs [concat $lDiseaseIDs [lsort -unique -increasing $lDs]]
    #########################################

    #######################################################
    puts "$CTDChemicalsDiseasesTSVFile is loading..."
    set F [open $CTDChemicalsDiseasesTSVFile]
    while {[gets $F Line]>=0} {
	if {[regexp {^\#} $Line]} {continue}
	
	set ChemicalID ""
	foreach {ChemicalName ChemicalID CasRN DiseaseName DiseaseID DirectEvidence} [lrange [split $Line "\t"] 0 5] {}
	if {[string equal $ChemicalID     ""]} {continue}
	if {[string equal $DirectEvidence ""]} {continue}
	if {[string equal $DirectEvidence "therapeutic"]} {continue}
	
	if {[info exists TabOBOID2OBOIDWithParents($DiseaseID)]} {
	    lappend lChemicalIDs    $ChemicalID
	    foreach DiseaseID $TabOBOID2OBOIDWithParents($DiseaseID) {
		set Tab($ChemicalID,$DiseaseID) 1
	    }
	}
    }
    close $F
    puts "$CTDChemicalsDiseasesTSVFile loaded!"
    set lChemicalIDs    [lsort -unique -increasing $lChemicalIDs   ]
    #######################################################
 
    #######################################################
    set OUTFile "${OUTDir}/CTD_chemical2CTD_diseases.txt"
    ::File::Open $OUTFile 100 "no"
    ::File::Save $OUTFile [join [concat [list "ID"] $lDiseaseIDs] "\t"]
    foreach ChemicalID $lChemicalIDs {
	set lNewLine [list $ChemicalID]
	foreach DiseaseID $lDiseaseIDs {
	    set Info 0
	    if {[info exists Tab($ChemicalID,$DiseaseID)]} {set Info 1}
	    lappend lNewLine $Info
	}
	::File::Save $OUTFile [join $lNewLine "\t"]
    }
    ::File::Close $OUTFile
    #######################################################

    return $OUTFile
}
proc CTD_genes_diseases {CTDGenesDiseasesTSVFile OBOCTDDiseasesFile HomoloGeneFile OUTDir} {
    #######################################################
    puts "$HomoloGeneFile is loading..."
    set F [open $HomoloGeneFile]
    while {[gets $F Line]>=0} {
	set GeneID ""
	foreach {HGID TAXID GeneID GeneSymb} [lrange [split $Line "\t"] 0 3] {}
	if {$GeneID == ""} {continue}

	set TabGeneID2HGID($GeneID) $HGID
    }
    close $F
    puts "$HomoloGeneFile loaded!"
    #######################################################

    #########################################
    set ListOfDiseaseIDs {}
    ::OBO::Load   $OBOCTDDiseasesFile
    foreach INIOBOID [ListOfNonRedundantElement [::OBO::Query $OBOCTDDiseasesFile Term]] {
	set TabOBOID2OBOIDWithParents($INIOBOID) {}
	foreach OBOID [ListOfNonRedundantElement [::OBO::Parents $OBOCTDDiseasesFile $INIOBOID 1]] {
	    set OBOName [lindex [OBO::Query $OBOCTDDiseasesFile Term $OBOID name] 0]
	    lappend TabOBOID2OBOIDWithParents($INIOBOID) "$OBOName (${OBOID})"
	    lappend ListOfDiseaseIDs "$OBOName (${OBOID})"
	}
    }
    ::OBO::Delete $OBOCTDDiseasesFile
    set ListOfDiseaseIDs [lsort -unique -increasing $ListOfDiseaseIDs]
    #########################################
    
    #######################################################
    set ListOfGeneIDs    {}
    puts "$CTDGenesDiseasesTSVFile is loading..."
    set F [open $CTDGenesDiseasesTSVFile]
    while {[gets $F Line]>=0} {
	if {[regexp {^\#} $Line]} {continue}
	
	set GeneID ""
	foreach {GeneSymbol GeneID DiseaseName DiseaseID DirectEvidence} [lrange [split $Line "\t"] 0 4] {}
	if {[string equal $GeneID         ""]} {continue}
	if {[string equal $DirectEvidence ""]} {continue}
	#if {[string equal $DirectEvidence "therapeutic"]} {continue}
	
	if {![info exists TabGeneID2HGID($GeneID)]} {continue}
	set GeneID $TabGeneID2HGID($GeneID)

	if {[info exists TabOBOID2OBOIDWithParents($DiseaseID)]} {
	    lappend ListOfGeneIDs    $GeneID
	    foreach DID $TabOBOID2OBOIDWithParents($DiseaseID) {
		set Tab($GeneID,$DID) 1
	    }
	}
    }
    close $F
    puts "$CTDGenesDiseasesTSVFile loaded!"
    set ListOfGeneIDs    [lsort -unique -increasing $ListOfGeneIDs   ]
    #######################################################

    #######################################################
    set OUTFile "${OUTDir}/CTD_HomoloGene2CTD_diseases.txt"
    ::File::Open $OUTFile 100 "no"
    ::File::Save $OUTFile [join [concat [list "ID"] $ListOfDiseaseIDs] "\t"]
    foreach GeneID $ListOfGeneIDs {
	set lNewLine [list $GeneID]
	foreach DiseaseID $ListOfDiseaseIDs {
	    set Info 0
	    if {[info exists Tab($GeneID,$DiseaseID)]} {set Info 1}
	    lappend lNewLine $Info
	}
	::File::Save $OUTFile [join $lNewLine "\t"]
    }
    ::File::Close $OUTFile
    #######################################################
    
    return $OUTFile
}

namespace eval ::OBO {
    variable TabOBOPanel
    proc Test {} {
	variable TabOBOPanel
	set OBOFile "/home/genouest/inserm625/fchalmel/GenomeAnnotation/GO/gene_ontology.obo"
	OBO::Load   $OBOFile
	foreach OBOID [OBO::Parents $OBOFile "GO:0051231" 1] {
	    puts "$OBOID -> [OBO::Query $OBOFile Term $OBOID name] -> [OBO::Query $OBOFile Term $OBOID namespace]"
	}
	foreach OBOID [OBO::Parents4GO $OBOFile [list "GO:0051231:IEA" "GO:0005524:EXP"] 1] {
	    set GOID [join [lrange [split $OBOID ":"] 0 1] ":"]
	    puts "$OBOID -> [OBO::Query $OBOFile Term $GOID name] -> [OBO::Query $OBOFile Term $GOID namespace]"
	}

	OBO::Delete $OBOFile

	return
    }
    proc Load {OBOFile} {
	variable TabOBOPanel
	if {![file exists $OBOFile]} {
	    puts stderr "$OBOFile does not exist!"
	    return 0
	}
	puts stdout "$OBOFile is loading..."
	
	set TabOBOPanel($OBOFile,LOADED) 1

	set InATermDefinition 0
	set InATypeDefinition 0
	set CurrentID         ""
	
	if {![info exists TabOBOPanel($OBOFile,Term)   ]} {set TabOBOPanel($OBOFile,Term)    {}}
	if {![info exists TabOBOPanel($OBOFile,Typedef)]} {set TabOBOPanel($OBOFile,Typedef) {}}
	
	set F [open $OBOFile]
	while {[gets $F Line]>=0} {
	    set Line   [string trim $Line]
	    set lLine  [split $Line ":"]
	    set HEADER    [lindex $lLine 0]
	    if {$HEADER == ""} {continue}
	    set REMAINDER [string trim [join [lrange $lLine 1 end] ":"]]
	    if {$HEADER == "\[Term\]"} {
		set InATermDefinition 1
		set InATypeDefinition 0
		set CurrentID         ""
		continue
	    }
	    if {$HEADER == "\[Typedef\]"} {
		set InATermDefinition 0
		set InATypeDefinition 1
		set CurrentID         ""
		continue
	    }
	    if {$InATermDefinition} {
		set TYPE "Term"
	    } elseif {$InATypeDefinition} {
		set TYPE "Typedef"
	    } else {
		if {![info exists TabOBOPanel($OBOFile,$HEADER)]} {set TabOBOPanel($OBOFile,$HEADER) {}}
		lappend TabOBOPanel($OBOFile,$HEADER) $REMAINDER
		continue
	    }
	    
	    if {$HEADER == "id"} {
		set     CurrentID                              [lindex [split $REMAINDER " "] 0]
		lappend TabOBOPanel($OBOFile,$TYPE)     $CurrentID
		set     TabOBOPanel($OBOFile,$TYPE,$CurrentID) $CurrentID
		continue
	    }
	    if {$CurrentID == ""} {continue}
	    
	    if {      $HEADER == "alt_id"      } {
		set ALTID [string trim [lindex [split $REMAINDER " "] 0]]
		set       TabOBOPanel($OBOFile,$TYPE,$ALTID) $CurrentID
		continue
	    } 
	    if {0 && $HEADER == "namespace"} {
		if {![info exists TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER) ]} {set TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER)  {}}
		lappend TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER)  [string trim $REMAINDER]
		continue
	    }
	    if {$HEADER == "is_a" || $HEADER == "relationship"} {
		set lREMAINDER [split $REMAINDER " "]
		if {$HEADER == "is_a"} {
		    set Relation  "is_a"
		    set RelatedID [string trim [lindex $lREMAINDER 0]]
		} else {
		    set Relation  [string trim [lindex $lREMAINDER 0]]
		    set RelatedID [string trim [lindex $lREMAINDER 1]]
		}
		if {$Relation == "is_a" || $Relation == "part_of"} {
		    if {![info exists TabOBOPanel($OBOFile,$TYPE,$CurrentID,Parents) ]} {set TabOBOPanel($OBOFile,$TYPE,$CurrentID,Parents)  {}}
		    if {![info exists TabOBOPanel($OBOFile,$TYPE,$RelatedID,Children)]} {set TabOBOPanel($OBOFile,$TYPE,$RelatedID,Children) {}}
		    lappend TabOBOPanel($OBOFile,$TYPE,$CurrentID,Parents)  $RelatedID
		    lappend TabOBOPanel($OBOFile,$TYPE,$RelatedID,Children) $CurrentID
		}
		if {![info exists TabOBOPanel($OBOFile,$TYPE,$CurrentID,Relationship)]} {set TabOBOPanel($OBOFile,$TYPE,$CurrentID,Relashionship) {}}
		lappend TabOBOPanel($OBOFile,$TYPE,$CurrentID,Relationship) [list $RelatedID $Relation]
		continue
	    }

	    if {![info exists TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER)]} {set TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER) {}}
	    lappend TabOBOPanel($OBOFile,$TYPE,$CurrentID,$HEADER) $REMAINDER
	    
	}
	close $F

	puts stdout "$OBOFile is loaded."
	return 1
    }
    proc Delete {{OBOFile ""}} {
	variable TabOBOPanel
	if {![info exists TabOBOPanel]} {return 0}

	if {$OBOFile != ""} {
	    foreach index [array names TabOBOPanel "$OBOFile,*"] {unset TabOBOPanel($index)}
	} else {
	    unset TabOBOPanel
	}
	return 1
    }
    proc Unset args {
	variable TabOBOPanel
	set ARG [join [lrange $args 0 end] ","]
	if {![info exists TabOBOPanel($ARG)]} {return 0}
	unset TabOBOPanel($ARG) 
	return 1
    }
    proc Query args {
	variable TabOBOPanel
	set ARG [join [lrange $args 0 end] ","]
	if {![info exists TabOBOPanel($ARG)]} {return}
	return [ListOfNonRedundantElement $TabOBOPanel($ARG)]
    }
    proc Exists args {
	variable TabOBOPanel
	set ARG [join [lrange $args 0 end] ","]
	if {![info exists TabOBOPanel($ARG)]} {return 0}
	return 1
    }
    proc Parents4GO {OBOFile ListOfOBOIDs WithChildren} {
	foreach OBOID $ListOfOBOIDs {
	    set lOBOID [split $OBOID ":"]
	    set GOID    [join [lrange $lOBOID 0 1] ":"]
	    set EviCode [lindex $lOBOID 2]
	    
	    if {$EviCode == ""} {set EviCode "NAS"}
	    if {![info exists Tab($EviCode)]} {set Tab($EviCode) {}}
	    lappend Tab($EviCode) $GOID
	}
	set ListOfParentOBOIDs {}
	foreach EviCode [array names Tab] {
	    foreach GOID [Parents $OBOFile $Tab($EviCode) $WithChildren] {
		set GOID "$GOID:${EviCode}"
		lappend ListOfParentOBOIDs $GOID
	    }
	}
	return $ListOfParentOBOIDs
    }
    proc Parents {OBOFile ListOfOBOIDs {WithChildren 1}} {
	set ListOfParentOBOIDs {}
	foreach OBOID [ListOfNonRedundantElement $ListOfOBOIDs] {
	    set OBOID [OBO::Query $OBOFile Term $OBOID]
	    foreach ParentOBOID [OBO::Parents_Recursif $OBOFile $OBOID] {
		if {[info exists TabDejaVu($ParentOBOID)]} {continue}
		set TabDejaVu($ParentOBOID) 1
		lappend ListOfParentOBOIDs $ParentOBOID
	    }
	}
	if {$WithChildren} {
	    return [ListOfNonRedundantElement [concat $ListOfOBOIDs $ListOfParentOBOIDs]]
	} else {
	    return [ListOfNonRedundantElement $ListOfParentOBOIDs]
	}
    }
    proc Parents_Recursif {OBOFile ParentOBOID {ListOfParentOBOIDs {}}} {
	foreach OBOID $ListOfParentOBOIDs {set TabDejaVu($OBOID) 1}
	foreach OBOID [OBO::Query $OBOFile Term $ParentOBOID Parents] {
	    set OBOID [OBO::Query $OBOFile Term $OBOID]
	    if {[info exists TabDejaVu($OBOID)]} {continue}
	    lappend ListOfParentOBOIDs $OBOID
	    set ListOfParentOBOIDs [OBO::Parents_Recursif $OBOFile $OBOID $ListOfParentOBOIDs]
	}
	return $ListOfParentOBOIDs
    }
    proc Children4GO {OBOFile ListOfOBOIDs WithParents} {
	foreach OBOID $ListOfOBOIDs {
	    set lOBOID [split $OBOID ":"]
	    set GOID    [join [lrange $lOBOID 0 1] ":"]
	    set EviCode [lindex $lOBOID 2]
	    
	    if {$EviCode == ""} {set EviCode "NAS"}
	    if {![info exists Tab($EviCode)]} {set Tab($EviCode) {}}
	    lappend Tab($EviCode) $GOID
	}
	set ListOfChildrenOBOIDs {}
	foreach EviCode [array names Tab] {
	    foreach GOID [Children $OBOFile $Tab($EviCode) $WithParents] {
		set GOID "$GOID:${EviCode}"
		lappend ListOfChildrenOBOIDs $GOID
	    }
	}
	return $ListOfChildrenOBOIDs
    }
 
    proc Children {OBOFile ListOfOBOIDs {WithParents 1}} {
	set ListOfChildOBOIDs {}
	foreach OBOID [ListOfNonRedundantElement $ListOfOBOIDs] {
	    set OBOID [OBO::Query $OBOFile Term $OBOID]
	    foreach ChildOBOID [OBO::Children_Recursif $OBOFile $OBOID] {
		if {[info exists TabDejaVu($ChildOBOID)]} {continue}
		set TabDejaVu($ChildOBOID) 1
		lappend ListOfChildOBOIDs $ChildOBOID
	    }
	}
	if {$WithParents} {
	    return [ListOfNonRedundantElement [concat $ListOfOBOIDs $ListOfChildOBOIDs]]
	} else {
	    return [ListOfNonRedundantElement $ListOfChildOBOIDs]
	}
    }
    proc Children_Recursif {OBOFile ParentOBOID {ListOfChildOBOIDs {}}} {

	foreach OBOID $ListOfChildOBOIDs {set TabDejaVu($OBOID) 1}
	foreach OBOID [OBO::Query $OBOFile Term $ParentOBOID Children] {
	    set OBOID [OBO::Query $OBOFile Term $OBOID]
	    if {[info exists TabDejaVu($OBOID)]} {continue}
	    lappend ListOfChildOBOIDs $OBOID
	    set ListOfChildOBOIDs [OBO::Children_Recursif $OBOFile $OBOID $ListOfChildOBOIDs]
	}
	return $ListOfChildOBOIDs
    }
}

namespace eval ::File {
    proc Open {File {N 100} {Complete "no"}} {
	variable Tab
	### if 'Complete' = yes, the information are added at the end of the file ###
	if {[string compare -nocase $Complete "yes"]} {set option "w"} else {set option "a"}
	set Channel          [open $File $option]
	set Tab($File)       $Channel
	set Tab($File,Lines) {}
	set Tab($File,N)     $N	
	return 1
    }
    proc Save {File Text {Complete "no"}} {
	### if 'Complete' = yes, the information are added at the end of the file ###
	variable Tab
	if {![info exists Tab($File)]} {
	    ::File::Open  $File 0     $Complete
	    ::File::Save  $File $Text $Complete
	    ::File::Close $File
	} else {
	    lappend Tab($File,Lines) $Text
	    if {$Tab($File,N) <= [llength $Tab($File,Lines)]} {
		puts $Tab($File)      [join $Tab($File,Lines) "\n"]
		set  Tab($File,Lines) {}
	    }
	}
	return $File
    }
    proc Close {File} {
	variable Tab
	if {[info exists Tab($File)] && [info exists Tab($File,Lines)]} {
	    if {0 < [llength $Tab($File,Lines)]} {
		puts $Tab($File) [join $Tab($File,Lines) "\n"]
	    }
	    close $Tab($File)

	    catch {unset Tab($File) Tab($File,Lines) Tab($File,N)}
	    return 1
	} else {
	    return 0
	}
    }
}
proc ListOfNonRedundantElement {Liste} {
   set ListOfNonRedundant {}
   foreach Elt $Liste {
       if {[info exists TabDejaVu($Elt)]} {continue}
       set     TabDejaVu($Elt)     1
       lappend ListOfNonRedundant $Elt
   }
   return $ListOfNonRedundant
}
CTD_genes_diseases     $CTDGenes2DiseasesFile     $CTDDiseaseOBOFile $HomoloGeneFile $OUTDir
CTD_chemicals_diseases $CTDChemicals2DiseasesFile $CTDDiseaseOBOFile $CTDChemicalsFile $EDDir $OUTDir
HomoloGene4R $HomoloGeneFile $OUTDir
exit