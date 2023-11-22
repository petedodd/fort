root <- plumber::plumb_api(package='fort', name='tbstatisticalserver')
plumber::pr_run(root, host="0.0.0.0", port=8080)