SELECT organism_id, common_name FROM organism
 WHERE EXISTS (
   SELECT * FROM feature
    WHERE feature.organism_id = organism.organism_id
      AND feature.type_id = (
      SELECT cvterm.cvterm_id
        FROM cvterm
        JOIN cv USING (cv_id)
       WHERE cv.name = 'sequence'
         AND cvterm.name = 'polypeptide'
      )
   )
