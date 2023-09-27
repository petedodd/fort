test_that("projection using scenario A returns I.mid of 290.58552", {
  # Arrange
  data <- jsonlite::fromJSON(system.file("plumber/tbstatisticalserver/projection.json",
                                         package = "fort"))
  set.seed(123)

  # Act
  
  result <- fort::projections(
    data$year,
    data$iHat,
    data$sEi,
    data$nHat,
    data$sEn,
    data$mHat,
    data$sEm,
    TXf = 0,
    HRd = 1,
    HRi = 1,
    ORt = 1
  )
  
  # Assert
  
  print(names(result$I.mid))
  print(result$I.mid[1])
  print(result$I.mid[2])
  
  testthat::expect_equal(result$I.mid[1], 290.58552)
})
