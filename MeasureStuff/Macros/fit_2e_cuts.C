

            if(vertexSec[0] < 5.7371 || vertexSec[0] > 5.8706 || vertexSec[1] < 5.7371 || vertexSec[1] > 5.8706)
            {
                continue;
            }
        ++ cut_counter[cc]; // cut 2 - vertex sector
        ++ cc;


        // Set the high and low energy index
        int highE_index = -1;
        int lowE_index = -1;
        if(electronEnergy[0] > electronEnergy[1])
        {
            highE_index = 0;
            lowE_index = 1;
        }
        else
        {
            highE_index = 1;
            lowE_index = 0;
        }



        // this was "old hotspots" changed to using new cuts based on
        // vertex position
        if(vertexInHotSpot[0])
        {
            hHotSpotLocation->Fill(vertexSec[0], vertexZ[0], 1.0);
        }
        if(vertexInHotSpot[1])
        {
            hHotSpotLocation->Fill(vertexSec[1], vertexZ[1], 1.0);
        }


                if( (electronBlockType[0] == 4) || (electronBlockType[1] == 4) ) continue; // not in petals
                if( (electronBlockType[0] == 5) || (electronBlockType[1] == 5) ) continue; // not in petals
        ++ cut_counter[cc]; // cut 3 - block type
        ++ cc;



                double cluster_total_energy = 0.;
                int cluster_hit_index = 0;
                int cesum = 0;

                //if(nGammaClusters > 0) continue;

                for(int i{0}; i < nGammaClusters; ++ i)
                {
                    /*
                    if(clusterEnergy[i] > 0.2)
                    {
                        ++ cesum;
                    }
                    */
                    /*
                    for(int j{0}; j < nInCluster[i]; ++ j)
                    {
                        if(clusterHitEnergy[cluster_hit_index] > 0.2)
                        {
                            ++ cesum;
                            cluster_total_energy += clusterHitEnergy[cluster_hit_index];
                        }
                        ++ cluster_hit_index;
                    }
                    */
                    //if(clusterEnergy[i] > 0.2) continue;
                    cluster_total_energy += clusterEnergy[i];
                }
                //if(cesum > 0) continue;
                //if(cluster_total_energy > 0.2) continue;
                if(cluster_total_energy > 0.0) continue;
                
        ++ cut_counter[cc]; // cut 4 - gamma energy
        ++ cc;




                if((trackSign[0] >= 0) || (trackSign[1] >= 0)) continue;
        ++ cut_counter[cc]; // cut 5 - track negative curvature
        ++ cc;
        

                if((eTrackLength[0] < 30. ) || (eTrackLength[1] < 30.)) continue;
        ++ cut_counter[cc]; // cut 6 - track length
        ++ cc;




                //if(vertexInHotSpot[0] || vertexInHotSpot[1]) continue;

                double vertexZ_1 = vertexZ[0];
                double vertexZ_2 = vertexZ[1];
                double vertexSec_1 = vertexSec[0];
                double vertexSec_2 = vertexSec[1];

                bool flag_continue = false;

                for(int elix = 0; elix < N_HOTSPOT_ELIPSE; ++ elix)
                {

                    double x1 = vertexSec_1;
                    double y1 = vertexZ_1;
                    double h1 = hotspot_elipse[elix][0];
                    double k1 = hotspot_elipse[elix][1];
                    double a1 = hotspot_elipse[elix][2];
                    double b1 = hotspot_elipse[elix][3];

                    double p1 = (std::pow(x1 - h1, 2.0) / pow(a1, 2.0)) +
                               (std::pow(y1 - k1, 2.0) / pow(b1, 2.0));

                    if(p1 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }

                    double x2 = vertexSec_2;
                    double y2 = vertexZ_2;
                    double h2 = hotspot_elipse[elix][0];
                    double k2 = hotspot_elipse[elix][1];
                    double a2 = hotspot_elipse[elix][2];
                    double b2 = hotspot_elipse[elix][3];

                    double p2 = (std::pow(x2 - h2, 2.0) / pow(a2, 2.0)) +
                               (std::pow(y2 - k2, 2.0) / pow(b2, 2.0));

                    if(p2 <= 1.0)
                    {
                        flag_continue = true;
                        break;
                    }
                               
                }

                if(flag_continue == true) continue;

        ++ cut_counter[cc]; // cut 7 (was cut 9) - track hot spot
        ++ cc;



                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    if(NAPromptGgHitsDist2Calo[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 0) continue;
        ++ cut_counter[cc];
        ++ cc;


                int ssum = 0;
                for(int i{0}; i < NAPromptGgHits; ++ i)
                {
                    if(NAPromptGgHitsDist2Vertex[i] < 15.)
                    {
                        ++ ssum;
                    }
                }
                if(ssum > 1) continue;
        ++ cut_counter[cc];
        ++ cc;
        

                int ssum = 0;
                if(electronSide[0] == electronSide[1])
                {
                    for(int i{0}; i < NAPromptGgHits; ++ i)
                    {
                        if(electronSide[0] != NAPromptGgHitsSide[i])
                        {
                            ++ ssum;
                        }
                    }
                }
                else if(electronSide[0] != electronSide[1])
                {
                    // tracks on opposite side
                    // do nothing in this case
                }
                if(ssum != 0) continue;
        ++ cut_counter[cc];
        ++ cc;


        
        
                if((firstGgHitLayer[0] > 0) || (firstGgHitLayer[1] > 0)) continue;
        ++ cut_counter[cc]; // cut 11 - Gg L0
        ++ cc;


                if(eeInternalProb < 0.01) continue;
                if(eeExternalProb > 0.01) continue;
        
        ++ cut_counter[cc]; // cut 12 (was cut 13) - prob
        ++ cc;
        

        double rPhi0 = (vertexR[0] * vertexSec[0] * TMath::TwoPi() / 20.);
        double rPhi1 = (vertexR[1] * vertexSec[1] * TMath::TwoPi() / 20.);
                if(fabs(vertexR[0] - vertexR[1]) > 4.) continue;
                if(fabs(vertexZ[0] - vertexZ[1]) > 8.) continue;
        ++ cut_counter[cc];
        ++ cc;
        

            if(mode_flag_2 == 1)
            {
                //if((electronEnergy[0] < 0.3) || (electronEnergy[1] < 0.3)) continue;
                if((electronEnergy[0] < 0.2) || (electronEnergy[1] < 0.2)) continue;
                // changed to 0.2 to permit application of systematics
            }
        ++ cut_counter[cc]; // cut 14 (was cut 15) - energy
        ++ cc;












            cut = true;
            if(thePhase == 0)
            {
                if((1869 <= run) && (run <= 3395)) cut = false;
            }
            else if(thePhase == 1)
            {
                if((3396 <= run) && (run <= 9186)) cut = false;
            }
            else if(nocutonphase == 1)
            {
                cut = false;
            }
            else
            {
                cut = true;
            }

            // accept P1 and P2 for testing purposes
            //cut = false;
            if(cut == true) continue;
        ++ cut_counter[cc]; // cut 1 - phase
        ++ cc;


