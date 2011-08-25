SELECT cvterm.name as cvterm_name, cv.name as cv_name
  FROM cvterm
  JOIN feature_cvterm on cvterm.cvterm_id = feature_cvterm.cvterm_id
  JOIN feature on feature.feature_id = feature_cvterm.feature_id
  JOIN cv on cv.cv_id = cvterm.cv_id
 WHERE (cv.name = 'RILEY' OR cv.name = 'genedb_products')
   AND feature.feature_id = %s
