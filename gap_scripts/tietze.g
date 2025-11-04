
# Example from GAP documentation: 48.6-2 SimplifyPresentation.
Test_0 := function()
    F2 := FreeGroup( "a", "b" );
    G := F2 / [ F2.1^9, F2.2^2, (F2.1*F2.2)^4, (F2.1^2*F2.2)^3 ];

    Print(G, " Size: ", Size(G), "\n");
    a := G.1;
    b := G.2;
    H := Subgroup( G, [ (a*b)^2, (a^-1*b)^2 ] );

    P := PresentationSubgroup( G, H );
    PrimaryGeneratorWords( P );
    TzOptions( P ).protected := 2; # Protects the first n generators from being removed. 
    TzOptions( P ).printLevel := 0; # From 0 to 3, less to more detail on print.

    TzPrintRelators( P );
    SimplifyPresentation( P );
    TzPrintRelators( P );
    return;
end;

# Example (a) from Havas, Kenne 1982 (??).
Test_1 := function()
    F2 := FreeGroup( "a", "b" );

    a := F2.1;
    b := F2.2;
    G := F2 / [ a^3, b^6, (a*b)^4, (a*b^2)^4, (a*b^3)^3, b^-2*a^-1*b^-2*a^-2*b^-2*a^-1*b^2*a*b^2*a^2*b^2*a];
    
    a := G.1;
    b := G.2;
    H := Subgroup( G, [ a, b^2 ] );

    Print("Index: ", Index(G, H), "\n");
    
    P := PresentationSubgroup( G, H );
    Print(P, "\n");
    
    #TzOptions( P ).protected := 2; # Protects the first n generators from being removed. 
    TzOptions( P ).printLevel := 0; # From 0 to 3, less to more detail on print.

    #TzPrintRelators( P );
    SimplifyPresentation( P );
    #TzPrintRelators( P );
    Print(P, "\n");
    return;
end;

# Example 4.8 from Linton 2025.
Test_2 := function()
    F2 := FreeGroup( "a", "b" );
    a := F2.1;
    b := F2.2;
    G := F2 / [ b^-1*a*b*a^-2, a^-1*b*a*b^-2 ];
    P := PresentationFpGroup(G);
    SimplifyPresentation( P );
    Print(P, "\n");
    return;
end;

Test_2();
