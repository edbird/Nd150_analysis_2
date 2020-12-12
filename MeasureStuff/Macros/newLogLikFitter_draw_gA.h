#ifndef NEWLOGLIKFITTER_DRAW_GA_H
#define NEWLOGLIKFITTER_DRAW_GA_H


bool is_even(const int var)
{
    if(var % 2 == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


double mean_x(const TH1D* const hist)
{
    double sum = 0.0;
    double psum = 0.0;
    for(Int_t i = 1; i <= hist->GetNbinsX(); ++ i)
    {
        // bin width = dx
        double width = hist->GetXaxis()->GetBinWidth(1);

        // probability p(x)
        double p = hist->GetBinContent(i);
        psum += p * width;

        // x
        double x = hist->GetXaxis()->GetBinCenter(i);

        sum += x * p * width;
        //sum += (x - mean) * p * width;
    }
//    std::cout << "psum=" << psum << std::endl;
    return sum;
}

double variance_x_fixed_mean(const TH1D* const hist, const double mean)
{
    double sum = 0.0;
    double psum = 0.0;
    for(Int_t i = 1; i <= hist->GetNbinsX(); ++ i)
    {
        // bin width = dx
        double width = hist->GetXaxis()->GetBinWidth(1);

        // probability p(x)
        double p = hist->GetBinContent(i);
        psum += p * width;

        // x
        double x = hist->GetXaxis()->GetBinCenter(i);

        // method 1
        sum += std::pow(x - mean, 2.0) * p * width;

        // method 2
        //sum += x * x * p * width;

        //sum += (x - mean) * p * width;
    }
//    std::cout << "psum=" << psum << std::endl;

    // method 2
    //sum -= mean * mean;

    return sum;
}

#endif // NEWLOGLIKFITTER_DRAW_GA_H
