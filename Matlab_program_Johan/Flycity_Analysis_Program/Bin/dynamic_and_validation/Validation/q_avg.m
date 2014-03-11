function [ q_avg ] = q_avg(q1, q2 ,q3, q4)


    % Function that caculates the average of a sequence of quaternions:
    
            M = length(q1);

            % Initial estimate:

            m = [q1(1); q2(1); q3(1); q4(1)];

            err = 1;

            while abs(err)<0.001

                err_t = ones(M,1);

                for k = 1:M

                    m_inv = [-m(1); -m(2); -m(3); m(4)]./(m(1)^2+m(2)^2+m(3)^2+m(4)^2);

                    m_inv = m_inv./norm(m_inv);

                    M_inv = [m_inv(4) m_inv(3) -m_inv(2) m_inv(1); ...
                             -m_inv(3) m_inv(4) m_inv(1) m_inv(2); ...
                             m_inv(2) -m_inv(1) m_inv(4) m_inv(3); ...
                             -m_inv(1) -m_inv(2) -m_inv(3) m_inv(4)];

                    M_q = M_inv*[q1(k); q2(k); q3(k); q4(k)];

                    M_q = M_q./norm(M_q);

                    err_t(k) = log(M_q);

                    clear m_inv M_inv M_q

                end

                err = sum(err_t);

                m = m*exp(err);

                m = m./norm(m);

                clear err_t

            end

            q_avg = [m(1); m(2); m(3); m(4)]/norm([m(1); m(2); m(3); m(4)]);    


end

