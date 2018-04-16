SELECT g.gene_id,e.left_position,e.right_position,g.strand 
FROM genes g INNER JOIN exons e USING(gene_id) 
INNER JOIN replicons r USING(replicon_id)
INNER JOIN genomes g USING(genome_id)
WHERE genome_id = 0
ORDER BY e.left_position ASC

/* 
 * Get all the genes within e coli sorted by their left position of the exon
 *
 */
SELECT g.gene_id,e.left_position,e.right_position,g.protein_id 
FROM genes g 
INNER JOIN exons e USING(gene_id) 
WHERE g.genome_id = 0
ORDER BY e.left_position ASC;


/**
 *
 * Get all the the orthologs to the protein id
 */
SELECT * FROM blast_1 b WHERE b.sseqid = 'WP_025591772' ORDER BY b.bitscore DESC;

/**
 *
 * Get gene information about the ortholog
 */
SELECT g.left_position, g.right_position, g.replicon_id from genes g WHERE g.protein_id = 'NP_414543'

/**
 * See how many genes exist between the orthologs
 */
SELECT * FROM genes g WHERE g.left_position > 583081 and g.right_position < 744629 and g.replicon_id = 2;
